% this code takes ~15 minutes on 2GHz processor

%% Reading_the_Metadata
% Get file lists
dirCohort1 = dir('Cohort1');
dirCohort2 = dir('Cohort2');

% Extract file names (including folders if present)
filesCohort1 = {dirCohort1.name}';
filesCohort2 = {dirCohort2.name}';

% Write to Excel
%writecell(filesCohort1, 'Cohort1_file_names.xlsx');
%writecell(filesCohort2, 'Cohort2_file_names.xlsx');

% Find overlapping file names
commonFiles = intersect(filesCohort1, filesCohort2);

%% Read the Meta Data Excel file
T = readtable('S4_Table_Chronic_Wound_Metadata.xlsx', ...
              'Sheet','Composite Data', ...
              'PreserveVariableNames', true);

%% 1) File names (column 1)
Metadata_filename = T{:,1};    % cell array of file names

%% 2) Week column: w1, w2, ... -> numeric
week_str = T.week;             % e.g. {'w1','w2',...}
week_num = str2double(erase(week_str,'w'));

%% 3) Patient column (cell array as is)
patient = T.patient;

%% 4) Healed column: cell + numeric label
healed = T.healed;             % cell array

healed_label = zeros(height(T),1);
healed_label(strcmp(healed,'Healer')) = 1;
healed_label(strcmp(healed,'Non-healer')) = 0;

%% 5) Treatment column: find variants + numeric labels
treatment = T.Treatment;       % cell array
[treatment_types,~,treatment_label] = unique(treatment);
% treatment_types : cell array of all variants
% treatment_label : numeric labels (1..N)

%% 6) Batch column: Batch1 / Batch2 -> numeric
batch = T.Batch;               % cell array
batch_label = zeros(height(T),1);
batch_label(strcmp(batch,'Batch1')) = 1;
batch_label(strcmp(batch,'Batch2')) = 2;

%% 7) "Healed by Visit #15": Yes / No -> numeric
healed_by_15 = T.("Healed by Visit #15");   % cell array
healed_by_15_label = zeros(height(T),1);
healed_by_15_label(strcmp(healed_by_15,'Yes')) = 1;
healed_by_15_label(strcmp(healed_by_15,'No'))  = 0;

% Indices where the labels differ
diff_idx = healed_label ~= healed_by_15_label;
diff_indices = find(diff_idx);

% % Analysis of differences between "healed" and "healed_by_15"
% samples_labelled_in_two_ways=patient(diff_indices);
% patients_labelled_in_two_ways=unique(patient(diff_indices));
% for i=1:numel(patients_labelled_in_two_ways)
%     n(i)=numel(find(ismember(patient, patients_labelled_in_two_ways(i))));
% end
% sum(n)-numel(samples_labelled_in_two_ways)
% k=find(ismember(patient,'BAART 095'));
% k1=healed_label(k);
% k2=healed_by_15_label(k);

%% 8) "Healed Visit #" column: read as cell array (mixed type)
healed_visit_number = T.("Healed Visit #");

healed_visit_label = nan(size(healed_visit_number));

for i = 1:numel(healed_visit_number)
    val = healed_visit_number{i};
    
    if ischar(val) || isstring(val)
        if strcmp(val,'Not healed')
            healed_visit_label(i) = NaN;
        else
            healed_visit_label(i) = str2double(val); % fallback
        end
    elseif isnumeric(val)
        healed_visit_label(i) = val;
    end
end
%% Compare meta file and actual filenames
% Merge two cohorts
mergedFiles = [filesCohort1; filesCohort2];

% Remove common elements
mergedFiles_noCommon = setdiff(mergedFiles, commonFiles, 'stable');
% Remove suffixes starting from _CKDL or _quant
mergedFiles_noCommon = regexprep( ...
    mergedFiles_noCommon, ...
    '(_CKDL.*|_quant.*)$', ...
    '' ...
);
overlap_with_metadata = intersect(mergedFiles_noCommon, Metadata_filename, 'stable');
only_in_cohorts = setdiff(mergedFiles_noCommon, Metadata_filename, 'stable');
only_in_metadata = setdiff(Metadata_filename, mergedFiles_noCommon, 'stable');

clear batch commonFiles diff_idx diff_indices dirCohort1 dirCohort2
clear filesCohort1 filesCohort2 healed healed_by_15 healed_visit_number
clear i k k1 k2 mergedFiles mergedFiles_noCommon n only_in_cohorts only_in_metadata
clear patients_labelled_in_two_ways samples_labelled_in_two_ways T
clear treatment val week_str overlap_with_metadata
%% Readig data
% Number of files
N = numel(Metadata_filename);

% Initialize outputs
Names_master = {};
NumReads_all = [];

for i = 1:N

    % Determine cohort folder
    if batch_label(i) == 1
        cohortFolder = 'Cohort1';
    elseif batch_label(i) == 2
        cohortFolder = 'Cohort2';
    else
        error('Unknown batch label at index %d', i);
    end

    % Get file prefix (e.g., 'D22')
    filePrefix = Metadata_filename{i};

    % Find matching file in the cohort folder
    fileStruct = dir(fullfile(cohortFolder, [filePrefix '_*.sf']));

    if isempty(fileStruct)
        error('No file found for prefix %s in %s', filePrefix, cohortFolder);
    elseif numel(fileStruct) > 1
        error('Multiple files found for prefix %s in %s', filePrefix, cohortFolder);
    end

    % Full file path
    filePath = fullfile(cohortFolder, fileStruct(1).name);

    % Read file
    T = readtable(filePath, 'FileType', 'text', 'Delimiter', '\t');

    % Extract Name and NumReads
    Names_current = T.Name;
    NumReads_current = T.NumReads;

    % First file: initialize master arrays
    if i == 1
        Names_master = Names_current;
        NumReads_all = NumReads_current;

    % Subsequent files: align by Name
    else
        % Preallocate new column with NaNs
        newNumReads = nan(numel(Names_master), 1);

        % Find matching indices
        [isMatch, idxInCurrent] = ismember(Names_master, Names_current);

        % Assign NumReads where names match
        newNumReads(isMatch) = NumReads_current(idxInCurrent(isMatch));

        % Append as new column
        NumReads_all(:, end+1) = newNumReads;
    end

end

%save(Readin_data_output.mat)
clear cohortFolder filePath filePrefix i idxInCurrent isMatch
clear N Names_current newNumReads NumReads_current T fileStruct
%% %Add gene names using transcripts
%load Readin_data_output.mat

GeneID=Names_master;

Tgl=readtable('converted_transcripts_to_genes.csv');
EnsID=Tgl{:,1};
GnNm=Tgl{:,3};

for i=1:numel(GeneID)
    tmp=GeneID(i);
    tmp1=string(tmp);
    tmp2=strfind(tmp1,'.');
    tmp3=char(tmp1);
    tmp4=tmp3(1:tmp2-1);
    tmp5(i,1)=cellstr(tmp4);
end
[tmp6,tmp7,tmp8]=intersect(EnsID,tmp5);
kk=zeros(size(GeneID));
for i=1:numel(tmp6)
    kk(tmp8(i))=tmp7(i);
end
clear tmp tmp1 tmp2 tmp3 tmp4 tmp5 tmp6 tmp7 tmp8

for i=1:numel(GeneID)
    if ~(kk(i)==0)
        Genes(i,1)=GnNm(kk(i));
    else
        Genes(i,1)={''};
    end
end
clear Tgl kk i GnNm EnsID

%remove rows with unknown genes:
idx = find(~cellfun(@isempty, Genes));
GeneID=GeneID(idx);
Genes=Genes(idx);
NumReads_all=NumReads_all(idx,:);
clear Names_master idx
%save('AddTranscriptsoutput.mat')

%% Use transcript with maximal expression as gene representative
%load AddTranscriptsoutput.mat

data=NumReads_all;
GenesU=unique(Genes);
for i=1:numel(GenesU)
    k=find(ismember(Genes,GenesU(i)));
    datK=log2(1+data(k,:));
        if numel(datK(:,1))>1
            tmp=sum(datK');
            [~,i1]=max(tmp);
        else
            i1=1;
        end
    dataU(i,:)=data(k(i1),:);
end

clear tmp i i1 k datK
clear data Genes GeneID NumReads_all
data=dataU;    clear dataU
Gene=GenesU; clear GenesU

%% Filter out genes with too many zeros
for i=1:numel(Gene)
    nz(i)=nnz(data(i,:));
end
k=find(nz>150);
data=data(k,:);
Gene=Gene(k);
clear i k nz
save ChronicWoundData.mat
