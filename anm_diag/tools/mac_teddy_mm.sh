#!/bin/bash
# teddy_mm 的 ANM 正式检验：所有必须在这台 Mac 上做的步骤，一条命令跑完。
#
#   bash <anm_diag 文件夹>/tools/mac_teddy_mm.sh
#
# 1. 下载 scPerturb 的四个文件（Zenodo 记录 10044268，约 1.7 GB，可断点续传），核对 Zenodo 上公布的 MD5。
# 2. 只看元数据（inspect，从不读计数矩阵），写到 runs/teddy_mm/inspect_<数据集>.txt。
# 3. 等 Claude 按元数据填好声明、过了检查，写入 runs/teddy_mm/GO（每 30 秒查一次）。
# 4. 每个数据集：封存 → 导出 token → TEDDY（MPS）→ 正式检验只跑一次，报告在 runs/<数据集>/anm_report.md。
#
# 中途断了，重跑同一条命令即可：做完的步骤会跳过；已封存的声明不会重封，已有结果的检验不会重跑。
# 进度和出错位置都记在 runs/teddy_mm/pipeline.log，当前阶段在 runs/teddy_mm/STATUS。
set -Eeuo pipefail      # -E: the failure trap also fires inside functions

if [ -z "${ANM_CAFFEINATED:-}" ] && command -v caffeinate >/dev/null 2>&1; then
  ANM_CAFFEINATED=1 exec caffeinate -i /bin/bash "$0" "$@"      # 跑的时候不让 Mac 闲置睡眠
fi

PKG="$(cd "$(dirname "$0")/.." && pwd)"
TEDDY_MM="${TEDDY_MM:-$(cd "$PKG/.." && pwd)}"   # anm_diag/ 所在的 teddy_mm 仓库；放在别处时用 TEDDY_MM=... 指定
PY="$TEDDY_MM/.venv/bin/python"
DATA="$TEDDY_MM/data/raw/perturb"
STATES="$TEDDY_MM/data/processed"
OUT="$PKG/runs/teddy_mm"
SUMS="$PKG/tools/scperturb_zenodo_10044268.md5"
URL="${ANM_SCPERTURB_URL:-https://zenodo.org/records/10044268/files}"   # 只有测试时才改
DATASETS="frangieh2021 papalexi2021"
# 重跑时给新的封存件和结果目录加后缀（比如 ANM_RUN_TAG=_headfix），旧的结果原样保留
TAG="${ANM_RUN_TAG:-}"
# 设计列（不是由数据算出来的列）：只看这些列的取值计数
SHOW="perturbation,condition,gene,guide_id,sgRNA,perturbation_type,cell_line,replicate"
export PYTHONUNBUFFERED=1

mkdir -p "$DATA" "$STATES" "$OUT"
LOG="$OUT/pipeline.log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }
status() { echo "$1" > "$OUT/STATUS"; }
trap 'log "失败：第 $LINENO 行：$BASH_COMMAND"; status FAILED' ERR

log "启动${TAG:+（重跑 ${TAG}）}：先找能用的 Python 环境，十几秒……"
# ---- 0. 找一个能用的 Python（要有 torch、anndata、scikit-learn、scipy、pyyaml、safetensors）
NEED="${ANM_PY_NEED:-import torch, anndata, sklearn, scipy, yaml, safetensors}"   # 只有测试时才改
SITE311="$TEDDY_MM/.venv/lib/python3.11/site-packages"
PYWAIT="${ANM_PY_WAIT:-90}"     # 每个候选最多等这么多秒（第一次加载 torch 可能要几十秒）；卡住就换下一个
# 用 perl 的 alarm 给每次尝试限时（macOS 没有 timeout 命令）；stdin 接 /dev/null，免得某个 Python 停下来等输入
try_py() { perl -e 'alarm shift; exec @ARGV' "$PYWAIT" "$1" -c "$NEED" < /dev/null > /dev/null 2>&1; }
try_py_site() { PYTHONPATH="$SITE311" perl -e 'alarm shift; exec @ARGV' "$PYWAIT" "$1" -c "$NEED" < /dev/null > /dev/null 2>&1; }
FOUND=""
if [ -n "${ANM_PYTHON:-}" ]; then
  log "  用指定的 Python：$ANM_PYTHON"
  if try_py "$ANM_PYTHON"; then FOUND="$ANM_PYTHON"
  elif try_py_site "$ANM_PYTHON"; then FOUND="$ANM_PYTHON"; export PYTHONPATH="$SITE311"
  else log "  指定的 Python 加载不了需要的包（或超过 ${PYWAIT} 秒）"; fi
fi
if [ -z "$FOUND" ]; then
  log "  试 teddy_mm 的 .venv"
  if try_py "$PY"; then FOUND="$PY"; fi
fi
if [ -z "$FOUND" ]; then
  # teddy_mm 的 .venv 链接到 Homebrew 的 python@3.11；那个解释器不在了，就找这台 Mac 上别的 Python 3.11，直接用 .venv 里装好的包
  for c in /opt/homebrew/opt/python@3.11/bin/python3.11 /opt/homebrew/Cellar/python@3.11/*/bin/python3.11 \
           /opt/homebrew/bin/python3.11 /usr/local/bin/python3.11 \
           /Library/Frameworks/Python.framework/Versions/3.11/bin/python3.11 \
           "$HOME"/.pyenv/versions/3.11*/bin/python3.11 "$HOME"/.local/share/uv/python/cpython-3.11*/bin/python3.11 \
           $(command -v python3.11 2>/dev/null || true); do
    [ -x "$c" ] || continue
    log "  试 $c"
    if try_py_site "$c"; then FOUND="$c"; export PYTHONPATH="$SITE311"; break; fi
  done
fi
PY="$FOUND"
if [ -z "$PY" ]; then
  # 最后一招：在 anm_diag 文件夹里另建一个环境（不动 teddy_mm 的 .venv），按 teddy_mm 的 requirements.txt 装包
  VENV="$PKG/.venv"
  if ! try_py "$VENV/bin/python"; then
    BASE=""
    for c in python3.13 python3.12 python3.11 python3.10 /opt/homebrew/bin/python3 /usr/local/bin/python3 python3; do
      p="$(command -v "$c" 2>/dev/null)" || continue
      if "$p" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)' >/dev/null 2>&1; then BASE="$p"; break; fi
    done
    if [ -z "$BASE" ]; then
      log "teddy_mm 的 .venv 用不了，这台 Mac 上也找不到 Python 3.10 以上的版本。先运行 brew install python@3.12，再重跑这条命令。"
      status FAILED; exit 1
    fi
    log "teddy_mm 的 .venv 用不了（它链接的 Homebrew python@3.11 不在了）。用 $BASE 在 $VENV 新建一个环境，按 teddy_mm 的 requirements.txt 装包，几分钟。"
    "$BASE" -m venv "$VENV"
    "$VENV/bin/python" -m pip install --quiet --upgrade pip 2>&1 | tee -a "$LOG"
    "$VENV/bin/python" -m pip install -r "$TEDDY_MM/requirements.txt" 2>&1 | tee -a "$LOG"
  fi
  try_py "$VENV/bin/python" || { log "新环境里还是缺包：$NEED"; status FAILED; exit 1; }
  PY="$VENV/bin/python"
fi
log "Python：$PY $("$PY" -c 'import sys; print(sys.version.split()[0])')"
cd "$PKG"
status RUNNING
log "开始（anm_diag 在 ${PKG}）"

md5_of() { if command -v md5 >/dev/null 2>&1; then md5 -q "$1"; else md5sum "$1" | cut -d' ' -f1; fi; }
fetch() { curl -L --fail --retry 5 --retry-delay 10 "$@" < /dev/null; }

# ---- 1. 下载并核对 MD5
status DOWNLOADING
while read -r want name; do
  [ -n "${name:-}" ] || continue
  f="$DATA/$name"
  if [ -f "$f" ] && [ "$(md5_of "$f")" = "$want" ]; then log "已下载：$name"; continue; fi
  log "下载：$name"
  if ! fetch -C - -o "$f" "$URL/$name?download=1"; then
    log "续传不成，整个重下：$name"
    [ -f "$f" ] && mv -f "$f" "$f.partial"
    fetch -o "$f" "$URL/$name?download=1"
  fi
  if [ "$(md5_of "$f")" != "$want" ]; then
    log "MD5 不符，整个重下一次：${name}（坏文件留作 $name.bad）"
    mv -f "$f" "$f.bad"
    fetch -o "$f" "$URL/$name?download=1"
    [ "$(md5_of "$f")" = "$want" ] || { log "MD5 仍然不符：$name"; status FAILED; exit 1; }
  fi
done < "$SUMS"
"$PY" tools/verify_downloads.py "$DATA" "$SUMS" | tee -a "$LOG"

# ---- 2. 只看元数据
status INSPECTING
for ds in $DATASETS; do
  f="$OUT/inspect_$ds.txt"
  if [ ! -s "$f" ]; then
    log "元数据：$ds"
    "$PY" -m anm_diag inspect --paths "paths/${ds}.yaml" --show "$SHOW" > "$f.tmp"
    mv "$f.tmp" "$f"
  fi
done

# ---- 3. 等声明填好
if [ ! -f "$OUT/GO" ]; then
  status WAITING_FOR_DECLARATIONS
  log "下载和元数据都好了。回到 Claude 的对话说一声“好了”：我按元数据填好、检查声明后写入 runs/teddy_mm/GO，这里会自己接着跑。别关这个窗口。"
  while [ ! -f "$OUT/GO" ]; do sleep 30; done
  log "收到 GO：$(cat "$OUT/GO")"
fi

# ---- 4. 先把所有声明一起封存（任何一个数据集出结果之前），再逐个 token → TEDDY → 正式检验（每个只跑一次）
for ds in $(cat "$OUT/GO"); do
  S="declarations/teddy_mm_${ds}${TAG}.sealed.yaml"; P="paths/${ds}.yaml"
  grep -q "${ds}_teddy_states.npz" "$P" || { log "路径文件 $P 的 teddy_states 不是 ${ds}_teddy_states.npz"; status FAILED; exit 1; }
  if [ ! -f "$S" ]; then
    status "SEALING $ds"; log "封存：$ds"
    "$PY" -m anm_diag seal "declarations/teddy_mm_${ds}.yaml" --paths "$P" --out "$S" 2>&1 | tee -a "$LOG"
  fi
done
for ds in $(cat "$OUT/GO"); do
  S="declarations/teddy_mm_${ds}${TAG}.sealed.yaml"
  P="paths/${ds}.yaml"
  TOK="$OUT/${ds}_tokens.npz"
  ST="$STATES/${ds}_teddy_states.npz"
  R="runs/${ds}${TAG}"
  if [ ! -f "$TOK" ]; then
    status "TOKENS $ds"; log "导出 token：$ds"
    "$PY" -m anm_diag export-tokens "$S" --adapter perturb_cite --paths "$P" --out "$TOK.tmp.npz" 2>&1 | tee -a "$LOG"
    mv "$TOK.tmp.npz" "$TOK"
  fi
  if [ ! -f "$ST" ]; then
    status "TEDDY $ds"; log "TEDDY 前向（只存 Z）：$ds"
    "$PY" tools/layer_states.py "$TOK" "$ST.tmp.npz" --z-only 2>&1 | tee -a "$LOG"
    mv "$ST.tmp.npz" "$ST"
  fi
  if [ ! -f "$R/anm_record.json" ]; then
    status "TESTING $ds"; log "正式检验（只跑一次）：$ds"
    "$PY" -m anm_diag run "$S" --adapter perturb_cite --paths "$P" --out "$R" 2>&1 | tee -a "$LOG"
  fi
  log "完成：${ds}，报告在 $PKG/$R/anm_report.md"
done
status DONE
log "全部完成。"
