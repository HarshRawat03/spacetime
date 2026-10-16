import sys, asyncio
from playwright.async_api import async_playwright
async def main(src, out, width=1280, scheme='light', height=900, full=True, wait=600):
    async with async_playwright() as p:
        b = await p.chromium.launch(executable_path='/opt/pw-browsers/chromium-1194/chrome-linux/chrome', args=['--no-sandbox'])
        pg = await b.new_page(viewport={'width':int(width),'height':int(height)}, color_scheme=scheme, device_scale_factor=1)
        logs=[]
        pg.on('console', lambda m: logs.append(f'{m.type}: {m.text}'))
        pg.on('pageerror', lambda e: logs.append(f'PAGEERROR: {e}'))
        await pg.goto(src)
        await pg.wait_for_timeout(int(wait))
        await pg.screenshot(path=out, full_page=(full=='1' or full is True))
        sw = await pg.evaluate('document.documentElement.scrollWidth')
        print('scrollWidth', sw, 'viewport', width)
        for l in logs: print(l)
        await b.close()
a=sys.argv
asyncio.run(main(a[1], a[2], a[3] if len(a)>3 else 1280, a[4] if len(a)>4 else 'light', a[5] if len(a)>5 else 900, a[6] if len(a)>6 else '1', a[7] if len(a)>7 else 600))
