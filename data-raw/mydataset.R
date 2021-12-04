## code to prepare `mydataset` dataset goes here

mydata = read.table("sample_data.txt", header = TRUE)

usethis::use_data(mydata, overwrite = TRUE)
