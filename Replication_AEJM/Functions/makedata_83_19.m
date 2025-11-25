% Upload the data from 1983Q1:2019Q4:

datafile = 'Data/US_bivariate.xlsx';
[num,txt] = xlsread(datafile);
dates_raw = txt(2:end,1);
startDate = '1983-01-01';
endDate = '2019-10-01';
start_ind = find(strcmp(startDate,dates_raw));
end_ind = find(strcmp(endDate,dates_raw));
data_excel = num(start_ind:end_ind,:);
raw_variables = txt(1,2:end);
dates = datetime(dates_raw(start_ind:end_ind));

for ii = 1:length(raw_variables)
    rawData.(raw_variables{ii}) = data_excel(:,ii);
end

%% Convert data

variable_list = {%  Short name, code generating variable, name,           transformation
                    'y',        'log(rawData.GDPC1)',     'GDP',          'diff';...
                    'p',        'log(rawData.GDPDEF)',    'Inflation',    'diff'};


variables = variable_list(:,1)'; % Short name of variable
var_legends = variable_list(:,3)'; % Descriptive name of variable

% Load data into the DATA object
for ii = 1:size(variable_list,1)
    codeString = strcat('DATA.',variable_list{ii,1},' = ',variable_list{ii,2},';');
    eval(codeString);
end

data    = nan(size(DATA.y,1)-1,size(variables,2)); % Data matrix

% Compute first differences (where specified)
for ii = 1:size(variables,2)
    if strcmp('diff',variable_list{ii,4})
        data(:,ii) = DATA.(variables{ii})(2:end,:)-DATA.(variables{ii})(1:end-1,:);
    else
        data(:,ii) = DATA.(variables{ii})(2:end,:);
    end
end

% Multiply by 100
data = data*100;

dates = dates(2:end);

[T,n] = size(data);