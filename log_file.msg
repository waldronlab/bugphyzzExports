Error in deparse(Exp[[i]]) != "" && is.call(curExp) : 'length = 5' in coercion to 'logical(1)' 
Traceback: 
 [1] "10: simplifyElementOp(curExp, useElementOp = FALSE, isTop = TRUE)"                                                                
 [2] "9: parserFunc(curLevel, codeMetaInfo, Exp)"                                                                                       
 [3] "8: ProcessCodeSingle(parserFunc, updateFunc, codePack$codeMetaInfo, "                                                             
 [4] "       curLevel, codePack$previousExp, codePack$processedExp, codePack$i, "                                                       
 [5] "       sub_Exp)"                                                                                                                  
 [6] "7: processCodePack(parserFunc, checkFunc, updateFunc, curLevel, "                                                                 
 [7] "       codePack, 0)"                                                                                                              
 [8] "6: parserFrame_hidden(parserFunc, checkFunc, updateFunc, codeMetaInfo)"                                                           
 [9] "5: parserFrame(RLevel1_parserFunc, RLevel1_checkFunc, RLevel1_updateFunc, "                                                       
[10] "       codeMetaInfo)"                                                                                                             
[11] "4: RParser1(codeMetaInfo0)"                                                                                                       
[12] "3: .compileGPUCode(FUN, parms, .macroParms = .macroParms, .options = .options)"                                                   
[13] "2: gpuSapply_singleDev(X, FUN, ..., .macroParms = .macroParms, .device = .device, "                                               
[14] "       .options = .options)"                                                                                                      
[15] "1: gpuSapply(X = data_ready, FUN = function(x) {"                                                                                 
[16] "       msg <- unique(x$Attribute_group)"                                                                                          
[17] "       msg <- paste0(\"Propagating \", msg, \"...\")"                                                                             
[18] "       log_print(msg)"                                                                                                            
[19] "       input_tbl <- complete(distinct(select(x, NCBI_ID, Attribute, "                                                             
[20] "           Score, Evidence)), NCBI_ID, Attribute, fill = list(Score = 0, "                                                        
[21] "           Evidence = \"\"))"                                                                                                     
[22] "       l <- split(input_tbl, factor(input_tbl$NCBI_ID))"                                                                          
[23] "       tree$Do(function(node) {"                                                                                                  
[24] "           if (!is.null(l[[node$name]])) {"                                                                                       
[25] "               node[[\"table\"]] <- l[[node$name]]"                                                                               
[26] "           }"                                                                                                                     
[27] "       })"                                                                                                                        
[28] "       tree$Do(asr, traversal = \"post-order\")"                                                                                  
[29] "       tree$Do(inh, traversal = \"pre-order\")"                                                                                   
[30] "       data_tree_tbl <- dplyr::filter(dplyr::relocate(dplyr::bind_rows(purrr::discard(tree$Get(function(node) node[[\"table\"]], "
[31] "           simplify = FALSE), ~all(is.na(.x)))), NCBI_ID), Evidence %in% "                                                        
[32] "           c(\"\", \"asr\", \"inh\") | is.na(Evidence))"                                                                          
[33] "       data_with_values <- bind_rows(data_tree_tbl, x)"                                                                           
[34] "       all_node_names <- tree$Get(function(node) node$name, simplify = TRUE)"                                                     
[35] "     ..."                                                                                                                         
Warning messages:
1: Missing columns in biosafety level. Missing columns are: Accession_ID
2: Missing columns in halophily. Missing columns are: Accession_ID
3: Missing columns in hemolysis. Missing columns are: Accession_ID
4: Missing columns in motility. Missing columns are: Accession_ID
5: Missing columns in pathogenicity human. Missing columns are: Accession_ID
6: Missing columns in spore formation. Missing columns are: Accession_ID
