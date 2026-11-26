const {chromium}=await import(process.env.PLAYWRIGHT_MODULE||'playwright');
const origin=process.env.CHALLENGE_QA_ORIGIN||'http://127.0.0.1:8792';
const browser=await chromium.launch({headless:true,executablePath:process.env.CHROMIUM_EXECUTABLE});
const page=await browser.newPage();const errors=[];page.on('pageerror',e=>errors.push(e.message));
for(const path of ['process/','review','activity','credits','research/']){
 for(const width of [390,1440]){
  await page.setViewportSize({width,height:1000});await page.goto(origin+'/challenge/'+path);await page.evaluate(()=>document.fonts.ready);
  if(await page.evaluate(()=>document.documentElement.scrollWidth>innerWidth))throw Error('Overflow '+path+' '+width);
  await page.screenshot({path:'/tmp/longevidence-'+path.replace('/','')+'-'+width+'.png',fullPage:true});
 }
 const links=await page.locator('a[href^="/challenge/"]').evaluateAll(a=>a.map(e=>e.getAttribute('href')));
 for(const link of new Set(links)){const r=await page.request.get(origin+link);if(!r.ok())throw Error('Broken '+link);}
}
await page.goto(origin+'/challenge/process/');await page.getByRole('link',{name:'Your credit',exact:true}).click();
if(!page.url().endsWith('#credit'))throw Error('Credit anchor failed');
if(errors.length)throw Error(errors.join('; '));
console.log('Process/review/activity/credits/research: desktop + mobile, anchors, local links and rendering passed.');
await browser.close();
