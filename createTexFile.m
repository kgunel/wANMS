function createTexFile(nVar,nBenchmarkFunc)

fileName = ['Results_' num2str(nVar)];
fid = fopen(['tex\' fileName '.tex'],'w+');

%fprintf(fid,'\\begin{table}[htbp]\n');
fprintf(fid,'\\centering{\n');
caption = [' \caption{The experimental results of NMS, ANMS and wANMS over 25 independent runs on $f_1 - f_{' num2str(nBenchmarkFunc)...
    '}$ benchmark functions of ' ...
    num2str(nVar) ' variables.}'];
label = ['\label{tbl:results_' num2str(nVar) '}'];
fprintf(fid,' \\begin{small}\n');
fprintf(fid,'  \\begin{longtable}{ccccc}\n');
fprintf(fid,'%s\n',[caption ' ' label ' \\']);
fprintf(fid,'     \\textbf{$f$-Id}  &  \\textbf{Results}  &  \\textbf{NMS}  &  \\textbf{ANMS}  & \\textbf{wANMS}  \\\\ \n');
fprintf(fid,'     \\hline\n');

for funId = 1:nBenchmarkFunc
  [funName, ~, ~, ~] = Cost(zeros(1,nVar),funId);
  dataFile = ['Result\' funName '_' num2str(nVar) '_Results.mat'];
  load(dataFile);
  
  [NMS_min_str, NMS_max_str, NMS_mean_str, NMS_std_str] = generateStats(NMS_Sol);
  [ANMS_min_str, ANMS_max_str, ANMS_mean_str, ANMS_std_str] = generateStats(ANMS_Sol);
  [wANMS_min_str, wANMS_max_str, wANMS_mean_str, wANMS_std_str] = generateStats(wANMS_Sol);
  
  fprintf(fid,'     \\multirow{4}{*}{$f_{%d}$} & $f_{\\text{Min}}$   & %s & %s & %s \\\\ \n',funId, NMS_min_str, ANMS_min_str, wANMS_min_str );
  fprintf(fid,'                              & $f_{\\text{Worst}}$ & %s & %s & %s \\\\ \n', NMS_max_str, ANMS_max_str, wANMS_max_str );
  fprintf(fid,'                              & $f_{\\text{Mean}}$  & %s & %s & %s \\\\ \n', NMS_mean_str, ANMS_mean_str, wANMS_mean_str );
  fprintf(fid,'                              & $f_{\\text{StDev}}$ & %s & %s & %s \\\\ \n', NMS_std_str, ANMS_std_str, wANMS_std_str );
  fprintf(fid,'     \\hline\n');
end
fprintf(fid,'   \\end{longtable}\n');
fprintf(fid,'  \\end{small}\n');
%fprintf(fid,'\\end{table}\n');
fprintf(fid,'}\n');
fclose(fid);
end