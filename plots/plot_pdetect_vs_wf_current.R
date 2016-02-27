
nbr_of_pt_spline = 1000

my_mse = function(v, w){
  norm(as.matrix(v-w))^2
}

dir_names = list.dirs(path = ".", full.names = FALSE, recursive = FALSE)

for(dir_name in dir_names){
  pdetect_file = paste(dir_name, paste("/", args[1], sep=""), sep="")
  weightfield_file = paste(dir_name, paste("/", args[2], sep=""), sep="")
  plot_file = paste(dir_name, paste("/", args[3], sep=""), sep="")

  pdetect_I_vs_t = read.table(pdetect_file, header=TRUE)
  wf_I_vs_t = read.table(weightfield_file, header=TRUE)

  time_col_name = "TIME.ns.."
  current_col_name = "Itot.uA.."

  #remove useless columns
  wf_I_vs_t = wf_I_vs_t[,c(time_col_name ,"Itot.uA..")]

  min_time = min(c(min(pdetect_I_vs_t[time_col_name ]), min(wf_I_vs_t[time_col_name ])))
  max_time = max(c(max(pdetect_I_vs_t[time_col_name ]), max(wf_I_vs_t[time_col_name ])))
  xlim = c(min_time, max_time)

  min_current = min(c(min(pdetect_I_vs_t[current_col_name]), min(wf_I_vs_t[current_col_name])))
  max_current = max(c(max(pdetect_I_vs_t[current_col_name]), max(wf_I_vs_t[current_col_name])))
  ylim = c(min_current, max_current)

  pdf(paste(plot_file,".pdf",sep=""))
  plot(wf_I_vs_t, col="red", xlim=xlim, ylim=ylim)
  points(pdetect_I_vs_t, col="green")
  legend(mean(xlim),mean(ylim),c("weightfield","pdetect"), lty=c(1,1),
  lwd=c(2.5,2.5), col=c("red","green"))
  dev.off()

  pdet_spline = splinefun(pdetect_I_vs_t)
  wf_spline = splinefun(wf_I_vs_t)

  t = seq(from=min_time, to=max_time, by=(max_time-min_time)/nbr_of_pt_spline)

  pdet_y = pdet_spline(t)
  wf_y = wf_spline(t)

  print(my_mse(pdet_y, wf_y))

  #Plot ratio
  ratio = pdet_y/wf_y
  pdf(paste(plot_file,"_ratio.pdf",sep=""))
  #print(ratio)
  plot(x=t,ratio, col="red")
  dev.off()
}
