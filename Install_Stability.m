
disp(' ');
disp('Building of combinatorial laplacian mex files...');

cd('Combinatorial');

mex stability_louvain_LCL.cpp community.cpp graph_binary.cpp;

cd ..;

disp('Building of normalised laplacian mex files...');

cd('Normalised');

mex stability_louvain_LNL.cpp community.cpp graph_binary.cpp;

cd ..;

disp('Moving build files to bin directory...');

cd('Combinatorial');
files=dir('stability_louvain_LCL.mex*');
for i=1:length(files)
    copyfile(files(i).name,'../bin');
end

cd('../Normalised');
files=dir('stability_louvain_LNL.mex*');
for i=1:length(files)
    copyfile(files(i).name,'../bin');
end

disp(' ');
disp('Do you wish to add stability to your Matlab path?');
disp('This will allow you to execute stability as a standard Matlab function.');
reply = input('Y/N [Y]: ', 's');


while ~strcmpi(reply(1),'Y') && ~strcmpi(reply(1),'N')
    reply = input('Please answer by ''Y'' or ''N''. Do you wish to add stability to your Matlab path? Y/N [Y]: ', 's');
end

if strcmpi(reply(1),'Y')
	cd('../bin');
	path(path,pwd);
    output_savepath = savepath;
    if ~output_savepath
        disp(' ');
		disp('Stability binaries sucessfully added to your Matlab path.');
    else
		warning('There was a problem adding stability to your Matlab path. You may need to run Matlab as a superuser. Alternatively you can save your path to a different location by calling SAVEPATH with an input argument that specifies the fullpath. For MATLAB to use that path in future sessions, save the path to ''pathdef.m'' in your MATLAB startup folder. ');
    end
end

cd ..;
disp(' ');
disp('Done');
