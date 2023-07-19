Error in file(file, "rt") : cannot open the connection to 'https://docs.google.com/spreadsheets/d/e/2PACX-1vQb6gDQdhCMWHLQdkzMM4bgTPPosYCZvjzJu4NyOYaSBGOCaOPoQh0OiUr4btnJ8EAWn7NiF-T9JaCE/pub?gid=0&single=true&output=csv' 
Traceback: 
[1] "6: file(file, \"rt\")"                                                 
[2] "5: read.table(file = file, header = header, sep = sep, quote = quote, "
[3] "       dec = dec, fill = fill, comment.char = comment.char, ...)"      
[4] "4: utils::read.csv(url)"                                               
[5] "3: dplyr::distinct(utils::read.csv(url))"                              
[6] "2: .importSpreadsheets(keyword = keyword)"                             
[7] "1: physiologies(phys_names, full_source = FALSE)"                      
Warning message:
In file(file, "rt") :
  URL 'https://docs.google.com/spreadsheets/d/e/2PACX-1vQb6gDQdhCMWHLQdkzMM4bgTPPosYCZvjzJu4NyOYaSBGOCaOPoQh0OiUr4btnJ8EAWn7NiF-T9JaCE/pub?gid=0&single=true&output=csv': status was 'SSL peer certificate or SSH remote key was not OK'
