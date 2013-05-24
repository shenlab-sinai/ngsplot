1) Copy "ngsplot", "multiplot" and "replot" folders to <your galaxy_install_dir>/tools/ directory.   

2) Edit  <your galaxy_install_dir>/tool_conf.xml and add the following lines in the <toolbox> section: 

   <section name="NGSPLOT Tools" id="ngs_plot">     
      <tool file="ngsplot/ngs.plot.xml"/>    
      <tool file="replot/replot.xml"/>     
      <tool file="multiplot/multiplot.xml"/>    
   </section> 

3) Restart Galaxy.
