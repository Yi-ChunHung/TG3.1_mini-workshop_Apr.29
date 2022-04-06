function wcal = ReadInput(inputfile)
    %# read lines
    fid = fopen(inputfile,'rt');
    C = textscan(fid, '%s', 'Delimiter',''); C = C{1};
    fclose(fid);
    
    %# start/end of each structure
    startIdx = find(ismember(C, 'GreenSur'));
    endIdx = find(ismember(C, 'endGreenSur'));
    
    %# array of strucutres
    N = numel(startIdx);
    wcal = struct('Ni',0, 'ef',0, 'Emu',0, 'Ne',0, 'np',0,'nn',0, ...
                  'kc',0, 'kp1',0, 'kp2',0, 'kp3',0, 'surface',0, 'hop_d',0, 'file','');
    wcal = repmat(wcal,[N 1]);
    
    %# parse and store each structure in the array
    for i=1:numel(startIdx)
        %# parse key/value of struct
        s = C(startIdx(i)+1:endIdx(i)-1); 
        s = regexp(s, '(\w+)\s*[:=]\s*([^%$]*)(?:%[^$]*)?', 'tokens', 'once');
        s = vertcat(s{:});
    
        %# try to parse as numbers
        v = str2double(s(:,2));
        s(~isnan(v),2) = num2cell(v(~isnan(v)));
    
        %# store: struct.key = value
        for j=1:size(s,1)
            wcal(i).(s{j,1}) = s{j,2};
        end
    end
    