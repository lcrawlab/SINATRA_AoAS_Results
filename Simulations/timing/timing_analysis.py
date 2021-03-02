import pandas as pd 

if __name__ == "__main__":
	timing_df = pd.read_csv('timings_sd.csv', index_col = 0)
	timing_df = timing_df.sort_values(['num.cones', 'dir.per.cone','ec.curve','num.shapes'], ascending=True)
	df = timing_df.set_index(['num.cones','dir.per.cone','ec.curve','num.shapes']) 
	print(df)


	with pd.option_context('display.max_rows', None, 'display.max_columns', None): 
    	print(df.unstack())
	
	

