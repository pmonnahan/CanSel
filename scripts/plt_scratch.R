vars = c("tajD", 'pi', 'Fst_Pat')
tdb = as.data.frame(tdb)

tdb = tbl(conn, '10_nstats')

tdb %<>% filter(pop %in% pops)

plotList = list()
for(i in 1:length(vars)){
    if(i==1){
       for(j in 1:22){
         tdb = as.data.frame(tbl(conn, paste0(j,'_nstats')))
         tdb %<>% filter(pop %in% pops)
         if(j==1){
            plt <- plotly_build(plot_ly(tdb, x = ~start, y = ~get(vars[i]), color = ~pop, legendgroup = ~pop, type = 'scatter', mode = 'lines', alpha = 0.5) %>% layout(yaxis = list(title = vars[i])))
         } else{
           plt %<>% add_trace(x = ~start, y = ~get(vars[i]), color = ~pop, legendgroup = ~pop, type = 'scatter', mode = 'lines', alpha = 0.5)         }
         }
    }
    else{
      for(j in 1:22){
        if(j==1){
          plt <- plotly_build(plot_ly(tdb, x = ~start, y = ~get(vars[i]), color = ~pop, legendgroup = ~pop, type = 'scatter', mode = 'lines', alpha = 0.5) %>% layout(yaxis = list(title = vars[i])))
        } else{
          plt %<>% add_trace(x = ~start, y = ~get(vars[i]), color = ~pop, legendgroup = ~pop, type = 'scatter', mode = 'lines', alpha = 0.5)            }
      }
    }
    plotList[[vars[i]]] = plt
  }
}

fig = subplot(plotList, nrows = length(plotList), shareX = TRUE)
# fig1 <- plot_ly(dat, x = ~start, y = ~tajD, color = ~pop, legendgroup = ~pop, type = 'scatter', mode = 'lines', alpha = 0.5)
# fig2 <- plot_ly(dat[dat$variable=='pi',], x = ~start, y = ~value, color = ~pop, legendgroup = ~pop, type = 'scatter', mode = 'lines', alpha = 0.5)
# fig <- subplot(figs, nrows=length(figs), shareX = TRUE)

button_list = list()
for (f in 1:22){
  new_list = list(method="restyle", args = list("visible", seq(1,22)==f), label = paste0("Chr", f))
  button_list = append(button_list, new_list)
}

fig <- fig %>% layout(
  xaxis = list(
    rangeslider = list(type = "numeric")),
  updatemenus = button_list))  

fig


# This is assuming that the visibility vector will be applied to all variables plotted.  May not be true.
bed = read.table("~/Documents/Research/ALL/ALL_genes.bed")
genebutton_list = list()
for (f in 1:nrow(bed)){
  gene = bed[f,]
  visible = seq(1,22)==gene[,1]) #Must make visible the chromosome track and the gene track.  Assuming gene tracks is are last n tracks on plot
  new_list = list(method="restyle", args = list(list("visible", visible, list(xaxis = list(range = c(gene[,2] - offset, gene[,3] + offset)))), label = gene[,4])
  genebutton_list = append(genebutton_list, new_list)
}

genefig <- fig %>% layout(updatemenus = genebutton_list))  

# Put all gene bars onto the top plot ie. vars[1]
height = max(tdb[vars[1]]) + (0.05 * max(tdb[vars[1]]))
for (f in 1:nrow(bed)){
  gene = bed[f,]
  pos = seq(gene[,2],gene[,3])
  plotList[[vars[1]]] <- add_trace(x = pos, y = rep(height, length(pos)), type = 'scatter', mode = 'lines', color = 'black')

updatemenus <- list(
  list(
    active = -1,
    type = 'buttons',
    buttons = list(
      list(
        label = '50s',
        method = "relayout",
        args = list(list(xaxis = list(range = c(gene[,2] - offset, gene[,3] + offset)))), 
      list(
        label = 'all',
        method = "relayout",
        args = list(list(xaxis = list(range = c(min(x), max(x))))))
    )
  )
)

p <- plot_ly(df, type = 'scatter', mode = 'markers') %>%
  add_trace(x = ~x, y = ~y1) %>%
  add_trace(x = ~x, y = ~y2)

p <- p %>% layout(
  title = "Button Restyle",
  xaxis = list(domain = c(0.1, 1)),
  yaxis = list(title = "y", range = c(0,10)),
  updatemenus = list(
    list(
      type = "buttons",
      y = 0.8,
      buttons = list(
        
        list(method = "restyle",
             args = list("visible", c(T,T)),
             label = "All"),
        
        list(method = "restyle",
             args = list("visible", c(F,T)),
             label = "Trace 0"),
        
        
        list(method = "restyle",
             args = list("visible", c(T,F)),
             label = "Trace 1")))
  ))



p <- p %>% layout(
  title = "Button Restyle",
  xaxis = list(domain = c(0.1, 1)),
  yaxis = list(title = "y", range = c(0,10)),
  updatemenus = list(
    list(
      type = "buttons",
      y = 0.8,
      buttons = list(
        list(method = "restyle",
             args = list("visible", c(T,T)),
             label = "Chr1"),
        
        list(method = "restyle",
             args = list("visible", c(F,T)),
             label = "Trace 0"),
        
        
        list(method = "restyle",
             args = list("visible", c(T,F)),
             label = "Trace 1")))
  ))