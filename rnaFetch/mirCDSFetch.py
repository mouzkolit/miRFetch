from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.action_chains import ActionChains
import time
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.support.select import Select
from rnaFetch.InitScrapper import InitScrapper
import pandas as pd
from selenium.common import exceptions
import warnings
from itertools import islice
import logging
from pySankey.sankey import sankey
from gprofiler import GProfiler
import matplotlib.pyplot as plt


class microTCDS(InitScrapper):
    """_summary_

    Args:
        InitScrapper (_type_): _description_
    """
    
    def __init__(self, search_dataframe, threshold = 0.9,  browser = "Chrome", headless = True):
        """_summary_

        Args:
            search_dataframe: can be list, dict, dataframe, set holoding the unique gene identifiers
            threshold (float, optional): _description_. Defaults to 0.9. should hold the threshold to select predictions
            browser (str, optional): _description_. Defaults to "Chrome". Determines which browser for fetching will be used, chrome will be downloaded so this is default
            headless (bool, optional): _description_. Defaults to True. Determines if browser is shown or not (headless False)

        Raises:
            ValueError: _description_
        """
        super().__init__(browser, headless)
        self.url = "https://dianalab.e-ce.uth.gr/html/dianauniverse/index.php?r=microT_CDS"
        self.driver.get(self.url)
        self.analysis = None
        self.searchable_dataframe = search_dataframe
        if not isinstance(threshold, float):
            raise ValueError("You provided the wrong data type as threshold")
        else:
            self._threshold = threshold
        self.data = None
        self.gprofiler = GProfiler(return_dataframe=True)
        
    def get_genes(self, gene_list):
        search_gene = " ".join(gene_list)
        return search_gene
    
    @property
    def threshold(self):
        """Getter Function of the setted threshold

        Returns:
            float: threshold height
        """

        print(f"Threshold is set to {self._threshold}")
        return self._threshold
    
    @threshold.setter
    def threshold(self, threshold):
        """Set the threshold to a float between 0-1

        Args:
            threshold (float): Threshold to be set before analysis
        """
        if self.analysis:
            if isinstance(threshold, float):
                if (threshold >= 0) and (threshold <= 1):
                    self._threshold = threshold
                    self.driver.find_element(By.ID, "adv_btn").click()
                    self.driver.find_element(By.ID, "threshold").clear()
                    self.driver.find_element(By.ID, "threshold").send_keys(str(threshold)) 
                    self.driver.find_element(By.NAME, "yt0").click()                      
                    print(f"Setted threshold successfully to {self._threshold}")
                else:
                    raise ValueError("Treshold not set between 0 and 1")  
            else:
                raise ValueError("Threshold is not of type float please specify correctly")    
        else:
            warnings.warn("Threshold can not selected before the Analysis! Analysis will be run at default 0.7 After the analysis you can change the threshold")

    def insert_search(self, genes):
        """Function to input the selected genes into the input box

        Args:
            genes (str): joined string of the imputed genes

        Returns:
            bool: True if run successfully otherwise connection error
        """
        input_area = self.driver.find_element(By.NAME, "keywords")
        input_area.send_keys(genes)
        input_area.send_keys(Keys.RETURN)
        time.sleep(0.5)
        check_state = self.remove_children()
        if check_state:
            progress_valid = self.check_progress()
            if progress_valid:
                return True
            else:
                raise ConnectionError("NO valid connection genes cannot be inserted into the Input Box")
    
    def remove_children(self):
        """Here we remove the not found genes using recursion

        Returns:
            bool: True if alle children were removed
        """
        if len(self.driver.find_elements(By.CLASS_NAME, "suggestions")) > 0:
            print("We found suggestions by microCDS, please check manually")
            get_element = self.driver.find_element(By.CLASS_NAME, "suggestions")
            get_children = get_element.find_elements(By.XPATH, "//form/p/input")
            try:
                for i in get_children:
                    i.click()
            except exceptions.StaleElementReferenceException as e:
                print("Deleted all suggestions")
            self.remove_children()
        else:
            return True
        
         
    def run_miRNA_analysis(self,threshold = 0.9, dictkey = "None"):
        """Start the miRNA analysis using the threshold and the searchable dataframe
        Args:
            threshold (float) <- threshold for the analysis between 0-1
            dictkey (str) <- if input data is of type dict give the list with the genes by the key in the dict
            
        returns:
            prediction_table(pd.DataFrame) <- DataFrame downloaded from microTCDS holding the predictions
        """
        prediction_table = pd.DataFrame()
     
        if not isinstance(dictkey, str):
            raise ValueError("You provided the wrong data type as dictkey please use a valid key as string")
        
        if not isinstance(threshold, float):
            raise ValueError("You provided wrong data type for threshold please use a float type")
        
        if isinstance(self.searchable_dataframe, pd.DataFrame):
            # This assumes a table that was provided by the biomart output and holds a column with the gene_ids
            gene_list = self.chunk(self.searchable_dataframe["Gene_ID"].astype(str).unique().tolist(), 200)
            
        if isinstance(self.searchable_dataframe, list):
            print("here we go")
            gene_list = {str(i) for i in self.searchable_dataframe}
            gene_list = self.chunk(gene_list, 200)
            
        if isinstance(self.searchable_dataframe, dict):
            gene_list = {str(i) for i in self.searchable_dataframe[dictkey]}
           
        # going through the individual chunks 
        for i in gene_list:
            genes_search = self.get_genes(list(i))
            self.insert_search(genes_search)
            self.analysis = 1
            self.threshold = threshold
            pred_table = self.load_prediction_table()
            prediction_table = pd.concat([prediction_table, pred_table])
            input_area = self.driver.find_element(By.NAME, "keywords")
            input_area.clear()
        self.data = prediction_table.copy()
        return prediction_table
        
            
    def check_progress(self):
        """Checks if the Job is still running or if all tables are already provided

         Returns:
            bool: if job was performed successful
        """
        print("Job is running")
        try:
            element = WebDriverWait(self.driver, 100).until(
                EC.presence_of_element_located((By.ID, "download_results_link")))
        except Exception as e:
            print(f"The following exception occured {e}")
            print("Check progress will be rerun!")
            self.check_progress()

        finally:
            print("All Data is queried successfully, Job completed")
            return True
            
    def load_prediction_table(self):
        """Download the table from the DIANA search page

        Args:
            element_to_check (Selenium Element): Check if the element of the download button exists

        Returns:
            pd.DataFrame: table with the prediction data from the table 
            
        """       
        print("Download of the Table initialized")
        download = self.driver.find_element(By.ID, "download_results_link")
        download_link = download.get_attribute("href")
        table = pd.read_csv(download_link, delimiter=",")
        table = table.dropna()
        table["Gene_ID"] = [i.split(" ")[0] for i in table["Gene Id(name)"]]
        table["Gene Symbol"] = [i.split(" ")[1].replace("(", "").replace(")","") for i in table["Gene Id(name)"]]
        
        print("Download Done and added to the Slot data_per_threshold")
        return table

        
    def get_mt_cds_overlap(self, microT_data, micro_CDS_data, groupby = "Sequence", number = 3):
        """Get the overlap between miRNA target space and queried RNA searches

        Args:
            microT_data (pd.DataFrame): Data received from the microT_data class
            micro_CDS_data (pd.DataFrame): Data fetched from the queried genes
            groupby (str, optional): _description_. Defaults to "Sequence".

        Returns:
            pd.DataFrame: merged_data is the whole merged table
            pd.DataFrame: grouped_table is the table grouped by Query and Mirna Name
        """
        merged_data = pd.merge(micro_CDS_data, microT_data, how = "left", right_on = "Gene_ID", left_on = "Gene_ID")
        grouped_table = pd.DataFrame(merged_data.groupby(["Query", "Mirna Name"])["Mirna Name"].count())
        grouped_table.columns = ["count"]
        grouped_table = pd.DataFrame(grouped_table.groupby(level = 0)["count"].nlargest(number)).reset_index(level = 1, drop = True)
        grouped_table = grouped_table.reset_index()
        grouped_table["count"] = grouped_table["count"]/20
        
        return merged_data, grouped_table
    
    def plot_sankey_mirna_overlap(self, sankey_dataframe):
        """Draw the Sankey Plot

        Args:
            sankey_dataframe (pd.DataFrame): Dataframe holding the data for the Sankey diagrom
        returns:
            fig (plt.Figure) <- holding the sankey plot
        """
             
        sankey(left=sankey_dataframe['Query'],
               right=sankey_dataframe['Mirna Name'],
               rightWeight=sankey_dataframe['count'],
               leftWeight = sankey_dataframe["count"],
               aspect=10,
               fontsize=10
               )
        # Get current figure
        fig = plt.gcf()

        # Set size in inches
        fig.set_size_inches(6, 12)

        # Set the color of the background to white
        fig.set_facecolor("w")
        
        return fig
        
    def get_gprofiles(self, data):
        """_summary_

        Args:
            data (_type_): _description_
        """
        print("needs to be implemetned")
        
        
    def chunk(self, it, size):
        """List will be chunked into equal size for request Since request size is only 1000

        Args:
            it (list): will be the iterable of the list
            size (int): describes the number of chunks based on size per chunk

        Returns:
            iter : iterable object for all chunks 
        """
        it = iter(it)
        return iter(lambda: tuple(islice(it, size)), ())
        

    

