% remove Intel directives from fortran codes

folder_name='compiled_on_Windows10c';

mkdir(folder_name);

homedir=pwd;

cd ..

MEXdir=[pwd '\MEX_files'];
cd(MEXdir);

filelist=dir('*.F')';

cd(homedir);

for file = filelist
    % remove Intel directives and write the fortran code to a new file
    FID=fopen([MEXdir '\' file.name]);
    C = fscanf(FID,'%c');
    fclose(FID);
    C=strrep(C,'!DEC$ simd','          ');
    FID=fopen(file.name,'w');
    fprintf(FID, '%c', C);
    fclose(FID);
    
    % compile the new file
    mex('-largeArrayDims', file.name, '-outdir',folder_name)
end


