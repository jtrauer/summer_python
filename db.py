from sqlalchemy import create_engine
import pandas as pd
import glob


class InputDB:
    """
    methods for loading input xls files
    """

    def __init__(self, database_name="Inputs.db"):
        """
        initialise sqlite database
        """
        self.database_name = database_name
        self.engine = create_engine("sqlite:///" + database_name, echo=False)

    def load_csv(self, input_path="xls/*.csv"):
        """
        load CSVs from inputPath
        """
        csv_file_list = glob.glob(input_path)
        for filename in csv_file_list:
            dataframe = pd.read_csv(filename)
            dataframe.to_sql(filename.split("\\")[1].split(".")[0], con=self.engine, if_exists="replace")

    def load_xslx(self, input_path ="xls/*.xlsx"):
        """
        load xslx from inputPath
        """
        excelFileList = glob.glob(input_path)
        available_sheets \
            = ["default_constants", "country_constants", "default_programs", "country_programs", "bcg_2014", "bcg_2015",
               "bcg_2016", "rate_birth_2014", "rate_birth_2015", "life_expectancy_2014", "life_expectancy_2015",
               "notifications_2014", "notifications_2015", "notifications_2016", "outcomes_2013", "outcomes_2015",
               "mdr_2014", "mdr_2015", "mdr_2016", "laboratories_2014", "laboratories_2015", "laboratories_2016",
               "strategy_2014", "strategy_2015", "strategy_2016", "diabetes", "gtb_2015", "gtb_2016", "latent_2016",
               "tb_hiv_2016", "spending_inputs"]

        for filename in excelFileList:
            xls = pd.ExcelFile(filename)

            if len(xls.sheet_names) == 1:
                df_name = xls.sheet_names[0]
                print(df_name)
            else:
                numSheets = 0
                while numSheets < len(xls.sheet_names):
                      sheet_name = xls.sheet_names[numSheets]
                      if sheet_name in available_sheets:
                          if sheet_name == "rate_birth_2015" or sheet_name == "life_expectancy_2015" :
                              df = pd.read_excel(filename, sheet_name=sheet_name, header = 3)
                          else:
                              df = pd.read_excel(filename, sheet_name=sheet_name)
                          print(sheet_name)
                          df.to_sql(sheet_name, con=self.engine, if_exists="replace")
                      numSheets = numSheets + 1

    def db_query(self, table_name, filter="", value="", column="*"):
        """
        method to query tablename
        """
        query = "Select " + column + " from  " + table_name
        if filter != '' and value != '':
            query = "Select " + column + " from  " + table_name + " Where " + filter + " = \'" + value + "\'"
        output_from_db = pd.read_sql_query(query, con=self.engine)
        return output_from_db


if __name__ == "__main__":

    input = InputDB()
    input.load_xslx()
    input.load_csv()
    res = input.db_query("bcg_2015", filter="Cname", value="Bhutan")
    print(res)
    # res = input.db_query("notifications_2016", filter="Country", value="Bhutan")
    # print(res)