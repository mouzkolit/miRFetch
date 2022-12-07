from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from webdriver_manager.chrome import ChromeDriverManager


class InitScrapper():
    def __init__(self, browser, headless = True):
        self.browser = browser
        self.driver = None
        self.options = Options()
        
        if headless:
            self.headless_state()
            

        if self.browser == "Chrome":
            self.chrome_init()
        elif self.browser == "Firefox":
            self.firefox_init()
        elif self.browser == "Edge":
            self.edge_init()

    def chrome_init(self):
        self.driver = webdriver.Chrome(ChromeDriverManager().install(), options = self.options)

    def firefox_init(self):
        self.driver = webdriver.Firefox()

    def edge_init(self):
        self.driver = webdriver.Edge()
        
    def headless_state(self):
        self.options.headless = True # Runs Chrome in headless mode.

