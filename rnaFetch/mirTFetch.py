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

class mirTFetch(InitScrapper):
    """
    Here we initialize the fetcher of the miRNA/tRNA prediction Data

    """

    def __init__(self, browser="Chrome", species: str = "Homo sapiens (Ensembl v84)") -> None:
        super().__init__(browser)

        self.url = "https://mrmicrot.imsi.athenarc.gr/?r=mrmicrot/index"
        self._threshold = 0.9
        self.driver.get(self.url)
        self.input_area = self.driver.find_element(By.NAME, "text_area")
        self.input_submit = self.driver.find_element(By.NAME, "yt0")
        self.species = species
        self._prediction_data = None
        self.utr_table = None

    def get_species_options(self):
        """Retrieve the Species that are availabe with the algorithm

        Returns:
            list: Selected Species listed which can be used for 
        """
        all_species_selected = []
        species = self.driver.find_element(By.NAME, "species")
        all_species = species.find_elements(By.XPATH, "//option")
        for option in all_species:
            all_species_selected.append(option.text)
        return all_species_selected

    def select_species(self, species):
        """_summary_

        Args:
            species (_type_): _description_
        """
        selected_species = Select(self.driver.find_element(By.NAME, "species"))
        selected_species.select_by_visible_text(species)
        print("Your species has been selected! You can now run the Pipeline with the following command \n self.run_miRNA_analysis(your_input)")
        time.sleep(0.8)

    @property
    def threshold(self):
        return self._threshold

    @threshold.setter
    def threshold(self, tres):
        """Setting the treshold for predictions

        Args:
            tres (float): Treshold
        """
        self._threshold = tres
        print(f"The current threshold is: {self._threshold}")

    @property
    def prediction_data(self):
        """Get prediction data from the analysis

        Returns:
            pd.DataFrame: Prediction Data
        """
        print("Retrieving the Prediction Data")
        if self._prediction_data.empty:
            print("No data were retrieved!")
        else:
            return self._prediction_data

    @prediction_data.setter
    def prediction_data(self, data):
        """Sets the prediction data into the slot

        Args:
            data (pd.DataFrame): Prediction Data Table
        """
        self._prediction_data = data

    def run_miRNA_analysis(self, nucleotide_dictionary):
        """The whole Analysis is performed by the mirT website
        Here the inputs to the websites will be mangaged and the
        state of the process defined

        Args:
            nucleotide_dictionary (_type_): _description_

        Returns:
            _type_: _description_
        """
        complete_table = pd.DataFrame()
        utr_full = pd.DataFrame()
        self.clear_input()
        for keys in nucleotide_dictionary:
            print(f"Started Submitting for the following key: {keys}")
            all_checks = []
            list_of_vales = nucleotide_dictionary[keys]
            analysis_started = self.insert_search(list_of_vales, keys)

            if analysis_started:
                self.input_submit.click()
                time.sleep(3)
                get_element = self.driver.find_element(
                    By.XPATH, "//div[@id='predictionResults']")
                all_children_by_xpath = get_element.find_elements(
                    By.XPATH, ".//div[@class='tarcloud_prediction']")

                for elements in all_children_by_xpath:
                    result = self.check_progress(elements)
                    all_checks.append(result)
                    
                print("All Data is queried successfully, Job completed")
                if all(all_checks) == True:
                    table, utr_table = self.load_prediction_table(
                        all_children_by_xpath)
                    table["Query"] = table.shape[0] * [keys]
                    utr_table["Query"] = utr_table.shape[0] * [keys]
                    complete_table = pd.concat([complete_table, table], axis=0)
                    utr_full = pd.concat([utr_full, utr_table], axis=0)
                    print("You can retrieve your table using self.prediction_data")

            else:
                print("Current Job was not succesfull")
                return None

            self.clear_input()

        self.utr_table = utr_full
        self.prediction_data = complete_table

    def insert_search(self, list_of_vales, keys):
        """Here the Search Terms will be inserted into the Text_Area
        of the Homepage

        Args:
            list_of_vales (list): list of Nucleotides

        Returns:
            bool: If Succesfully True will be returned
        """

        for values in list_of_vales:
            values = values.upper()
            check_status = self.check_nucleotides(values)

            if check_status:
                self.input_area.send_keys(values)
                ActionChains(self.driver).key_down(Keys.SHIFT).key_down(
                    Keys.ENTER).key_up(Keys.SHIFT).key_up(Keys.ENTER).perform()

            else:
                print(
                    "You didn't provide Nucleotides Sequence; should only contain (A,T,G,C) letters")
                return None
        print(f"Succesfully checked nucleotides for submission of key : {keys}")
        return True

    def load_prediction_table(self, element_to_check):
        """_summary_

        Args:
            element_to_check (_type_): _description_

        Returns:
            _type_: _description_
        """
        final_table = pd.DataFrame()
        utr_table_full = pd.DataFrame()
        for i in element_to_check:
            download = i.find_element(By.ID, "mrmicrot_result_download_button")
            download_link = download.get_attribute("href")
            table = pd.read_csv(download_link, delimiter="|", header=None)
            utrs, utr_table = self.get_utrs(download_link)
            table = table.dropna()
            table.columns = ["Sequence", "Transcript_ID", "Score"]
            table["utrs"] = utrs
            table = table[table["Score"] >= self.threshold]
            final_table = pd.concat([final_table, table], axis=0)
            utr_table_full = pd.concat([utr_table_full, utr_table], axis=0)
        return final_table, utr_table_full

    def get_utrs(self, download_link):
        data = pd.read_csv(download_link, header=None)
        list_of_utr_counts = []
        utr_count = None
        gene_list = []
        for i in data.iterrows():
            if utr_count == None:
                utr_count = 0
                current_gene = i[1][0].split("|")[1]
                continue
            if ">" in i[1][0]:
                current_gene = i[1][0].split("|")[1]
                list_of_utr_counts.append(utr_count)
                utr_count = 0
            else:
                utrs = i[1][0].split("\t")
                utrs.append(current_gene)
                gene_list.append(utrs)
                utr_count += 1

        list_of_utr_counts.append(utr_count)

        utr_dataframe = pd.DataFrame(
            gene_list, columns=["UTR", "Transcript Position", "Score", "Transcript_ID"])
        return list_of_utr_counts, utr_dataframe

    def check_progress(self, element_to_check):
        """Checks if the Job is still running or if all tables are already provided

        Args:
            element_to_check (_type_): the html element which we wait for

        Returns:
            bool: if job was performed successful
        """
        print("Job is running")
        try:
            element = WebDriverWait(element_to_check, 100).until(
                EC.presence_of_element_located((By.ID, "mrmicrot_result_download_button")))
        except Exception as e:
            print(f"The following exception occured {e}")
            print("Check progress will be rerun!")
            self.check_progress(element_to_check)

        finally:
            return True

    def check_nucleotides(self, nucleotides):
        """_summary_

        Args:
            nucleotides (str): str of nucleotides
        """
        nucleotides_set = set(nucleotides)
        if nucleotides_set.issubset({"A", "T", "C", "G"}):
            return True
        else:
            print(
                f"Here are the nucleotides you provided: {nucleotides_set}, please check, not supported by the Server!")

    def clear_input(self):
        """_summary_
        """
        self.input_area.clear()
        print("Input Area of Server cleared after Job")
