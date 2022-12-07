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

    def __init__(self, gene_list) -> None:
        self.chunked_list = gene_list
        self.workers = 2 * multiprocessing.cpu_count() + 1
        self.biomart_data = None
        self.merged_dataframe = None
        self.configure_data
        # here at a certain point more filters and attributes can be added

    def configure_data(self):
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
        new_list = self.chunk(self.chunked_list, 500)
        checked_list = self.progress_biomart(new_list)
        for i in checked_list:
            print(len(i["external_gene_name"]), len(
                i["Gene_ID"]), len(i["Transcript_ID"]))
            data = pd.DataFrame(i)
            final_annotation_dataframe = pd.concat(
                [final_annotation_dataframe, data], axis=0)
        self.biomart_data = final_annotation_dataframe.copy()
        return final_annotation_dataframe

    def progress_biomart(self, new_list):
        """Runs the Worker Thread for each chunk and returns an iterator of the list
        results

        Args:
            new_list (_type_): _description_

        Returns:
            _type_: _description_
        """
        with ThreadPoolExecutor(max_workers=self.workers) as executor:
            return executor.map(self.check_biomart, new_list, timeout=60)

    def check_biomart(self, chunked_list):
        """_summary_

        Args:
            chunked_list (_type_): _description_

        Returns:
            _type_: _description_
        """
        gene_dict = {"external_gene_name": [],
                     "Gene_ID": [], "Transcript_ID": []}
        chunked_list = list(chunked_list)
        server = "http://rest.ensembl.org"
        ext = "/lookup/id"
        search_gene = '{ "ids" :' + " " + f"{json.dumps(chunked_list)}" + " }"
        headers = {"Content-Type": "application/json",
                   "Accept": "application/json"}
        r = requests.post(server+ext, headers=headers, data=search_gene)
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        decoded = r.json()
        for keys in decoded:
            try:
                gene_dict["external_gene_name"].append(
                    decoded[keys]["display_name"][:-4])
                gene_dict["Gene_ID"].append(decoded[keys]["Parent"])
                gene_dict["Transcript_ID"].append(decoded[keys]["id"])
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

