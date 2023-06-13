c <- read.csv(file="casos-notifi.csv")

c$Month <- lubridate::ymd(c$Month)

  ggplot(c)+
    geom_col(aes(Month, Total_muestras), fill="#7AD151FF")+
    geom_col(aes(Month, Casos), fill="#440154FF") +
    geom_point(aes(Month, brotes*20), size=0.4, color="#2a788eff")+
    scale_y_continuous(
    sec.axis = sec_axis(trans = ~. /20, name="Number of outbreaks")) + 
    labs(fill="",
       x=NULL,
       y="Montly processed samples")+
    theme_minimal()+
    theme(text = element_text(size = 14))


  
  # Python by chatGTP
  
  import pandas as pd
  import matplotlib.pyplot as plt
  import seaborn as sns
  from matplotlib.ticker import FuncFormatter
  
  # Load the data
  c = pd.read_csv('casos-notifi.csv')
  
  # Convert Month to datetime format
  c['Month'] = pd.to_datetime(c['Month'])
  
  # Plot the data
  plt.figure(figsize=(10, 6))
  
  # Total_muestras
  plt.bar(c['Month'], c['Total_muestras'], color='#7AD151FF')
  
  # Casos
  plt.bar(c['Month'], c['Casos'], color='#440154FF')
  
  # brotes
  plt.scatter(c['Month'], c['brotes'] * 20, s=4, color='#2a788eff')
  
  # Secondary axis for number of outbreaks
  def outbreaks_formatter(x, pos):
    return int(x / 20)
  
  ax = plt.gca()
  ax2 = ax.twinx()
  ax2.yaxis.set_major_formatter(FuncFormatter(outbreaks_formatter))
  ax2.set_ylabel('Number of outbreaks')
  
  # Set labels and title
  plt.xlabel(None)
  plt.ylabel('Monthly processed samples')
  plt.title('Monthly Cases')
  
  # Adjust layout and display the plot
  plt.tight_layout()
  plt.show()
  