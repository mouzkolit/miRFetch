from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from webdriver_manager.chrome import ChromeDriverManager


class InitScrapper():
    
    def __init__(self, browser: str, headless: bool = True):
        """_summary_

        Args:
            browser (str): The name of the browser you would like to use e.g. Chrome/Firefox/Edge
            headless (bool, optional): If browser window should be shown or not. Defaults to True.
        """
        self.browser = browser
        self.driver = None
        self.options = Options()
        self.headless = headless
        
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
        self.assertEqual
        self.driver = webdriver.Firefox()

    def edge_init(self):
        self.driver = webdriver.Edge()
        
    def headless_state(self):
        self.options.headless = True # Runs Chrome in headless mode.

