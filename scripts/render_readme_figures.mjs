// Render the README figures (docs/assets/readme/*.png) from docs/assets/readme/src/figures.html.
// Needs Node 18+ and Playwright with a Chrome/Chromium: `npm i -D playwright` (or set PLAYWRIGHT_MODULE
// to a playwright install path). Usage: node scripts/render_readme_figures.mjs
import path from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const src = pathToFileURL(path.join(root, 'docs/assets/readme/src/figures.html')).href;
const outDir = path.join(root, 'docs/assets/readme');

const mod = await import(process.env.PLAYWRIGHT_MODULE ? pathToFileURL(process.env.PLAYWRIGHT_MODULE).href : 'playwright');
const { chromium } = mod.default ?? mod;
const browser = await chromium.launch({ channel: process.env.PW_CHANNEL || 'chrome' });

for (const theme of ['light', 'dark']) {
  const page = await browser.newPage({ viewport: { width: 1100, height: 900 }, deviceScaleFactor: 2 });
  await page.goto(`${src}?theme=${theme}`, { waitUntil: 'networkidle' });
  await page.evaluate(() => document.fonts.ready);
  await page.addStyleTag({ content: 'body{background:transparent!important}' });
  for (const id of ['hero', 'why', 'how']) {
    await page.locator(`#${id}`).screenshot({ path: path.join(outDir, `${id}-${theme}.png`), omitBackground: true });
    console.log('wrote', `${id}-${theme}.png`);
  }
  await page.close();
}
await browser.close();
