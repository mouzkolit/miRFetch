import unittest
from rnaFetch.mirCDSFetch import microTCDS


class MyTestCase(unittest.TestCase):

    @classmethod  
    def setUpClass(cls): 
        cls.micro = microTCDS(["hsa-miR-21-5p"])
        cls.micro_2 = microTCDS(1)

    @classmethod
    def tearDownClass(cls):
        """ Tears down the driver for completness"""
        cls.micro.driver.close()
        cls.micro_2.driver.close()

    def test_headless(self):
        """ Tests the state of the current initialized browser"""
        self.assertEqual(self.micro.headless, True, "Wrong Default Parameter")

    def test_input_data_type(self):
        """ Tests the input data to the microTCDS class"""
        self.assertRaises(TypeError, lambda: self.micro_2.run_miRNA_analysis())
        self.assertRaises(ValueError, lambda: microTCDS(["hsa-miR-21-5p"], threshold = "hello"))

    def test_list_mirnas_inserted(self):
        """ Check that there are no side effects after running the function to the base searchable
        dataframe"""
        self.micro.run_miRNA_analysis()
        self.assertEqual(self.micro.searchable_dataframe, ["hsa-miR-21-5p"])
        self.assertIsInstance(self.micro.searchable_dataframe, list)

    def test_output(self):
        """ Tests the output of the analysis 
        Should be a pd.DataFrame with 108 rows and 6 columns"""
        table = self.micro.run_miRNA_analysis()
        self.assertEqual(table.shape[0], 108)
        self.assertEqual(table.shape[1], 6)

if __name__ == '__main__':
    unittest.main()