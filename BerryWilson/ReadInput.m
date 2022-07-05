function wcal = ReadInput(inputfile)
    %# read lines
    fid = fopen(inputfile,'rt');
    C = textscan(fid, '%s', 'Delimiter',''); C = C{1};
    fclose(fid);
    
    %# start/end of each structure
    startIdx = find(ismember(C, 'WilsonLoop'));
    endIdx = find(ismember(C, 'endWilsonLoop'));
    
    %# array of strucutres
    N = numel(startIdx);
    wcal = struct('nVB',0, 'ib1',0, 'ib2',0, 'ib3',0, 'mode',0, 'ref','','ncpu',0);
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
    