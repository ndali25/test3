%% Load data

startDate = '01.01.2001';
endDate  = '01.03.2023';


[num,txt] = xlsread('Data/EuroArea_FRED.xlsx');

dates = txt(2:end,1);

startInd = find(strcmp(startDate,dates))-12;
endInd = find(strcmp(endDate,dates));

rawData.HICP = num(startInd:endInd,2);
rawData.INPR = num(startInd:endInd,1);
dates = dates(startInd+12:endInd);



%% Convert data

variable_list = {% Short name, code generating variable, name, transformation
    'y',    'log(rawData.INPR)',                            'GDP',              'diff';...
    'p',    'log(rawData.HICP)',                            'Prices',           'yoy'};


variables = variable_list(:,1)'; % Short name of variable
var_legends = variable_list(:,3)'; % Descriptive name of variable

% Load data into the DATA object
for ii = 1:size(variable_list,1)
    codeString = strcat('DATA.',variable_list{ii,1},' = ',variable_list{ii,2},';');
    eval(codeString);
end

data    = nan(size(DATA.y,1)-12,size(variables,2)); % Data matrix
% dates = x2mdate(rawData.GDPC1.Data(2:end,1)-693960,0,'datetime');


% Take first differences (where specified)
for ii = 1:size(variables,2)
    if strcmp('diff',variable_list{ii,4})
        data(:,ii) = DATA.(variables{ii})(13:end,:)-DATA.(variables{ii})(12:end-1,:);
    elseif strcmp('yoy',variable_list{ii,4})
        data(:,ii) = DATA.(variables{ii})(13:end,:)-DATA.(variables{ii})(1:end-12,:);
    else
        data(:,ii) = DATA.(variables{ii})(13:end,:);
    end
end

data = data*100;

[T,n] = size(data);