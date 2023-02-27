import pandas as pd
import requests
import sys
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
from itertools import islice
import multiprocessing
import warnings


class GeneBiomart():

    def __init__(self, gene_list, species = "homo_sapies", symbol = False) -> None:
        self.chunked_list = gene_list
        self.workers = 2 * multiprocessing.cpu_count() + 1
        self.biomart_data = None
        self.merged_dataframe = None
        self.configure_data
        self.species = species
        self._symbol = symbol
        # here at a certain point more filters and attributes can be added

    @property
    def symbol(self):
        return self._symbol
    
    def configure_data(self):
        """Identifies the number of workers retrieving the biomart from the api
        """
        print(f"The Biomart Tool was configures with: Nr of Workers: {self.workers}")
        
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

    def query_data(self):
        """Runs through the input list provide chunks of the list and spawns worker which 
        should submit list to the Biomart server

        Returns:
            _type_: _description_
        """
        final_annotation_dataframe = pd.DataFrame()
<<<<<<< HEAD
        chunked_gene_list = self.chunk(self.gene_list, 500)
        checked_list = self.progress_biomart(chunked_gene_list)

        for result in checked_list:
            data = pd.DataFrame(result)
            final_annotation_dataframe = pd.concat([final_annotation_dataframe, data], axis=0)
=======
        new_list = self.chunk(self.chunked_list, 100)
        checked_list = self.progress_biomart(new_list)
        for i in checked_list:
            print(len(i["external_gene_name"]), len(
                i["Gene_ID"]), len(i["Transcript_ID"]))
            data = pd.DataFrame(i)
            final_annotation_dataframe = pd.concat(
                [final_annotation_dataframe, data], axis=0)
>>>>>>> cb1144a793ed4c86001b0067776ac226667338eb
        self.biomart_data = final_annotation_dataframe.copy()
        return final_annotation_dataframe

    def progress_biomart(self, new_list: list):
        """Runs the Worker Thread for each chunk and returns an iterator of the list
        results

        Args:
            new_list (list): chunked list that is the input for the biomart retrieval

        Returns:
            _type_: _description_
        """
        with ThreadPoolExecutor(max_workers=self.workers) as executor:
            return executor.map(self.check_biomart, new_list, timeout=60)

    def check_biomart(self, chunked_list:list):
        """Retrieves the data from biomart either from symbol or from 

        Args:
            chunked_list (_type_): _description_

        Returns:
            _type_: _description_
        """
        
        chunked_list = list(chunked_list)
        
        if self.symbol:
            server = "https://rest.ensembl.org"
            ext = f"/lookup/symbol/{self.species}"
            search_gene = '{ "symbols" :' + " " + f"{json.dumps(chunked_list)}" + " }"
            print(search_gene)
        else:
            server = "http://rest.ensembl.org"
            ext = "/lookup/id"
            search_gene = '{ "ids" :' + " " + f"{json.dumps(chunked_list)}" + " }"
            
        headers = {"Content-Type": "application/json",
                   "Accept": "application/json"}
        
        try:
            r = requests.post(server+ext, headers=headers, data=search_gene)
            r.raise_for_status()
        except requests.exceptions.HTTPError as err:
            print(err)
            return gene_dict
        
        decoded = r.json()
        if self.symbol:
            return self.create_dict_from_decoding_symbol(decoded)
        else:
            return self.create_dict_from_decoding_ids(decoded)
        
    def create_dict_from_decoding_ids(self, decoded):
        """_summary_

        Args:
            decoded (_type_): _description_

        Returns:
            _type_: _description_
        """
        gene_dict = {"external_gene_name": [],
                     "Gene_ID": [], "Transcript_ID": []}
        for keys in decoded:
            try:
                gene_dict["external_gene_name"].append(
                    decoded[keys]["display_name"][:-4])
                gene_dict["Gene_ID"].append(decoded[keys]["Parent"])
                gene_dict["Transcript_ID"].append(decoded[keys]["id"])
            except:
                pass
        return gene_dict

    def create_dict_from_decoding_symbol(self, decoded):
        """_summary_

        Args:
            decoded (_type_): _description_

        Returns:
            _type_: _description_
        """
        gene_dict = {"external_gene_name": [],
                     "Gene_ID": [], "Transcript_ID": []}
        for keys in decoded:
            try:
                gene_dict["external_gene_name"].append(
                    keys)
                gene_dict["Gene_ID"].append(decoded[keys]["id"])
                gene_dict["Transcript_ID"].append(decoded[keys]["canonical_transcript"])
            except:
                pass

        return gene_dict

    def merge_annotation(self, data_input, data_constructed):
        """_summary_

        Args:
            data_input (_type_): _description_
            data_constructed (_type_): _description_
        """
        final_dataframe = pd.merge(
            data_input, data_constructed, how="left", left_on="Transcript_ID", right_on="Transcript_ID")
        self.merged_dataframe = final_dataframe.copy()
        return final_dataframe

