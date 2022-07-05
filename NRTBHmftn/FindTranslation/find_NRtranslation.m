% Codes that help to find the translation vectors of nano-ribbon
clear all

load slab_info_1.mat
hkl_1 = hkl;
BR1   = BR;
laod slab_info_2.mat
hkl_2 = hkl;
BR2   = BR;

%% the bases of hkl plane
t    = cell(1,4);
t{1} = BR1(1,:);
t{2} = BR1(2,:);
t{3} = BR2(1,:);
t{4} = BR2(2,:);

%% find the direction that perpendicular to both planes
for i = 1:4
    p1 = t{i}*hkl_1';
    p2 = t{i}*hkl_2';
    if abs(p1)<1e-4 && abs(p2)<1e-4
        j = i; % the desired direction 
        break
    end
end

%% find the translation vectors for narrow ribbon
p  = zeros(1,2);
if j>2
    j2 = setdiff([3,4],j); % translation vector on the hkl_2
    % find the j1: translation vector on the hkl_1
    for i = 1:2
        p(i) = t{i}*t{j}';
    end
    [~,j1] = min(p);
else
    j1 = setdiff([1,2],j); % translation vector on the hkl_1
    % find the j2: translation vector on the hkl_2
    for i = 1:2
        p(i) = t{i+2}*t{j}';
    end
    [~,j2_temp] = min(p);
    j2 = j2_temp + 2;
end

fprintf('The translation vectors: \n')
t1 = t{1}
t2 = t{2}
