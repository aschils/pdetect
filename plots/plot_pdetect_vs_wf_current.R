
nbr_of_pt_spline = 1000

my_mse = function(v, w){
  diff = v-w
  sum(diff*diff)
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

  min_time_pdetect = min(pdetect_I_vs_t[time_col_name ])
  min_time_wf = min(wf_I_vs_t[time_col_name ])
  min_time = min(c(min_time_pdetect, min_time_wf))
  max_time_pdetect = max(pdetect_I_vs_t[time_col_name ])
  max_time_wf = max(wf_I_vs_t[time_col_name ])
  max_time = max(c(max_time_pdetect, max_time_wf))
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

  min_time_ratio = max(c(min_time_pdetect, min_time_wf))
  max_time_ratio = min(c(max_time_pdetect,max_time_wf))

  t = seq(from=min_time_ratio, to=max_time_ratio, by=(max_time_ratio-min_time_ratio)/nbr_of_pt_spline)

  pdet_y = pdet_spline(t)
  wf_y = wf_spline(t)

  print(dir_name)
  #print(my_mse(pdet_y, wf_y))

  #Plot ratio
  ratio = (pdet_y-1)/(wf_y-1)
  pdf(paste(plot_file,"_ratio.pdf",sep=""))
  #print(ratio)
  plot(x=t,ratio, col="red", xlab = "Time (ns)")
  dev.off()
}
