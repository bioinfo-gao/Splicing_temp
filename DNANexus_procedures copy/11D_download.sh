# the link will expire after 24h. if not finished, have to download the remaining files
# https://linuxize.com/post/how-to-run-linux-commands-in-background/  
 wget -i sample_url.txt 
 wget -i 2022_07_26_remaining_url.txt 
 nohup wget -i 2022_07_26_remaining_url1.txt & # 2 3 4 5