import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data = pd.read_csv('/camhpc/ngs/projects/TST12086/dnanexus/20230305034844_Zhen.Gao/DSG_ranking_and_join/countries.csv')
data.head()

data_2007 = data[data.year == 2007]
datasort = data_2007.sort_values('population', ascending = False)
datasort = datasort.head(10)
datasort
x = range(10)
#The below code will create two plots. The parameters that .subplot take are (row, column, no. of plots).
plt.subplot(2,1,1)
#This will create the bar graph for poulation
pop = plt.bar(x, datasort['population']/10**6)
plt.ylabel('Population in Millions')
plt.xticks([],[])
#The below code will create the second plot.
plt.subplot(2,1,2)
#This will create the bar graph for gdp i.e gdppercapita divided by population.
gdp =plt.bar(x, datasort['gdpPerCapita'] * datasort['population'] / 10 ** 9)
plt.ylabel('GDP in Billions')
plt.xticks(x, datasort['country'], rotation='vertical')
plt.show()
