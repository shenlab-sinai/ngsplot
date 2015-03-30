1) Copy "ngsplot_intro", "ngsplot_main" and "ngsplot_replot" folders to <your galaxy_install_dir>/tools/ directory.   

2) Edit  <your galaxy_install_dir>/config/tool_conf.xml and add the following lines in the <toolbox> section: 

   <section name="NGSPLOT Tools" id="ngs_plot">     
      <tool file="ngsplot_intro/intro.xml"/>    
      <tool file="ngsplot_main/ngsplot.xml"/>    
      <tool file="ngsplot_replot/replot.xml"/>     
   </section> 

3) Restart Galaxy.
