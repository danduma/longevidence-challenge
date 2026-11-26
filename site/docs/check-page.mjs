const { chromium } = await import(process.env.PLAYWRIGHT_MODULE || 'playwright');
const origin = process.env.CHALLENGE_QA_ORIGIN || 'http://127.0.0.1:8767';
const output = process.env.CHALLENGE_QA_OUTPUT || '/tmp';
const browser = await chromium.launch({headless:true,executablePath:process.env.CHROMIUM_EXECUTABLE});
const context = await browser.newContext({permissions:['clipboard-read','clipboard-write']});
const page = await context.newPage();
const errors = [];
page.on('pageerror', e => errors.push(e.message));
page.on('console', m => { if(m.type() === 'error') errors.push(m.text()); });
await page.setViewportSize({width:1440,height:1000});
await page.goto(origin+'/challenge/');
await page.evaluate(async () => { await document.fonts.ready; await Promise.all(document.getAnimations().filter(a => a.effect?.getTiming().iterations !== Infinity).map(a => a.finished)); });
await page.screenshot({path:output+'/desktop.png',fullPage:true});
for (const width of [320,390,620,768,1024,1440]) {
  await page.setViewportSize({width,height:900});
  await page.evaluate(async () => { await Promise.all(document.getAnimations().filter(a => a.effect?.getTiming().iterations !== Infinity).map(a => a.finished)); });
  if(await page.evaluate(() => document.documentElement.scrollWidth > innerWidth)) throw Error('Overflow '+width);
  if(width === 390) await page.screenshot({path:output+'/mobile.png',fullPage:true});
}
await page.locator('#copy-prompt').click();
const copied = await page.evaluate(() => navigator.clipboard.readText());
if(!copied.includes(origin+'/challenge/join.txt')) throw Error('Wrong instruction URL');
await page.locator('summary').filter({hasText:'Where does the data go?'}).click();
if(!await page.locator('details[open]').count()) throw Error('FAQ did not open');
const links = await page.locator('a[href^="/challenge/"]').evaluateAll(es => es.map(e => e.getAttribute('href')));
for(const url of new Set(links)) {
  const r = await page.request.get(origin+url);
  if(!r.ok()) throw Error('Broken '+url);
}
const task = await (await page.request.get(origin+'/challenge/tasks/lc-001.json')).json();
if(task.id !== 'LC-001') throw Error('Missing task');
await page.evaluate(() => Object.defineProperty(navigator,'clipboard',{value:{writeText:async() => { throw Error('clipboard disabled'); }}}));
await page.locator('#copy-prompt').click();
if(!(await page.locator('#copy-status').textContent()).includes('manually')) throw Error('Fallback failed');
if(errors.length) throw Error(errors.join('; '));
console.log('Six viewport widths, desktop/mobile screenshots, fonts, copy/fallback, FAQ, instruction links and task checked.');
await browser.close();
