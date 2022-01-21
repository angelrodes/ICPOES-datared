%% Clear everything
clear
clc
close all hidden

%% emulate NORMPDF function locally to avoid using the package "statistics" in octave
normpdf_local = @(x,mu,sigma) exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);

%% Define Input parameters

% first we need to convert Qtegra data to standard ICP input

% input: csv with semicolon as delimiter
% ;;Raw.Average;(...);Raw.STD;(...)(;ExtCal.StandardConcentration
% ;;0;0;0;...
% ;;K 766.490 {44} (Radial);K 769.896 {44} (Radial);Ca 393.366 {86} (Radial);...
% 1;5% HNO3 MMR;-7.8025000000005091;...


% output: csv with the following columns:
%   Sample,Location,Analyte,Element,Wavelength,
%   ICPsignal,ICPsignal_error,NOMINAL_CONCENTRATIONS

% headers we are interested in:
xstring='ExtCal.StandardConcentration';
ystring='Raw.Average';
dystring='Raw.STD';
dystring2='Raw.RSD'; % just in case you did not include the previous one
delimiter=';';


% select file
xls_files = dir(fullfile('*.csv'));
XLS_files = dir(fullfile('*.CSV'));
xls_str = fliplr([{xls_files.name}]);
[s,v] = listdlg('Name','Input data','PromptString','Select a .csv file:',...
    'SelectionMode','single',...
    'ListString',xls_str);
profilefile=xls_str(s);
profilefile=profilefile{1};
outputfilename=[profilefile(1:end-4) '_INPUT.txt']; % this is the name of the new csv file we are going to generate

% read file
fileName=profilefile;
fid = fopen(fileName,'r');   %# Open the file
lineArray = cell(100,1);     %# Preallocate a cell array (ideally slightly larger than needed)
lineIndex = 1;               %# Index of cell to place the next line in
nextLine = fgetl(fid);       %# Read the first line from the file
while ~isequal(nextLine,-1)         %# Loop while not at the end of the file
    lineArray{lineIndex} = nextLine;  %# Add the line to the cell array
    lineIndex = lineIndex+1;          %# Increment the line index
    nextLine = fgetl(fid);            %# Read the next line from the file
end
fclose(fid);                 %# Close the file
lineArray = lineArray(1:lineIndex-1);  %# Remove empty cells, if needed
for iLine = 1:lineIndex-1              %# Loop over lines
    lineData = textscan(lineArray{iLine},'%s',...  %# Read strings
        'Delimiter',delimiter);
    lineData = lineData{1};              %# Remove cell encapsulation
    if strcmp(lineArray{iLine}(end),delimiter)  %# Account for when the line
        lineData{end+1} = '';                     %#   ends with a delimiter
    end
    lineArray(iLine,1:numel(lineData)) = lineData;  %# Overwrite line data
end
data=lineArray; % all the data in the csv file

% check if we have columns
if size(data,2)<4
    error(['Only ' num2str(size(data,2)) ' columns found: check that delimiter is "' delimiter '"'])
end



% check where is xstring
if sum(strcmp(ystring,data(1,:)))>1
    % correct position
elseif sum(strcmp(ystring,data(:,1)))>1
    % transpose
    data=data';
else
    error(['No ' xstring ' found in first row or column'])
end

% get data locations,samplenames,analytes,elements,wavelengths
sampleindex=5:size(data,1); % start looking for row headers on row 5
signalcolumns=find(strcmp(ystring,data(1,:)));
analytes=data(3,signalcolumns);


%% Build new database & save XXXX-INPUT.csv file
%   Sample,Location,Analyte,Element,Wavelength,
%   ICPsignal,ICPsignal_error,NOMINAL_CONCENTRATIONS
delete(outputfilename)
diary(outputfilename) % record what we are showing with "disp" and place in the new csv file
disp(['Sample,Location,Analyte,Element,Wavelength,ICPsignal,ICPsignal_error,NOMINAL_CONCENTRATIONS','RSD']) % print headers
linecount=0;
for row=sampleindex
    analytecount=0;
    for Analyte=analytes
        linecount=linecount+1;
        analytecount=analytecount+1;
        C = strsplit(Analyte{:},' ');
        element=C(1);
        wavelength=C(2);
        signalcolumn=find(strcmp('Raw.Average',data(1,:))&strcmp(Analyte,data(3,:)));
        dsignalcolumn=find(strcmp('Raw.STD',data(1,:))&strcmp(Analyte,data(3,:)));
        dsignalcolumn2=find(strcmp('Raw.RSD',data(1,:))&strcmp(Analyte,data(3,:)));
        nominalcolumn=find(strcmp('ExtCal.StandardConcentration',data(1,:))&strcmp(Analyte,data(3,:)));
        stringoutput=[data{row,2} ',' data{row,1} ',' Analyte{:} ',' element{:} ',' wavelength{:} ','...
            data{row,signalcolumn} ',' data{row,dsignalcolumn} ',' data{row,nominalcolumn} ',' data{row,dsignalcolumn2}];
        disp(stringoutput) % print line
    end
end
diary off % stop writing the csv file

%% Open XXXX-INPUT.csv file
% load file
% [num, txt, raw]= xlsread(profilefile);
% data=read_mixed_csv(profilefile,',');
fileName=outputfilename;
delimiter=',';
fid = fopen(fileName,'r');   %# Open the file
lineArray = cell(100,1);     %# Preallocate a cell array (ideally slightly
%#   larger than is needed)
lineIndex = 1;               %# Index of cell to place the next line in
nextLine = fgetl(fid);       %# Read the first line from the file
while ~isequal(nextLine,-1)         %# Loop while not at the end of the file
    lineArray{lineIndex} = nextLine;  %# Add the line to the cell array
    lineIndex = lineIndex+1;          %# Increment the line index
    nextLine = fgetl(fid);            %# Read the next line from the file
end
fclose(fid);                 %# Close the file
lineArray = lineArray(1:lineIndex-1);  %# Remove empty cells, if needed
for iLine = 1:lineIndex-1              %# Loop over lines
    lineData = textscan(lineArray{iLine},'%s',...  %# Read strings
        'Delimiter',delimiter);
    lineData = lineData{1};              %# Remove cell encapsulation
    if strcmp(lineArray{iLine}(end),delimiter)  %# Account for when the line
        lineData{end+1} = '';                     %#   ends with a delimiter
    end
    lineArray(iLine,1:numel(lineData)) = lineData;  %# Overwrite line data
end
data=lineArray; % our data is now in the right format! Bye bye Qtegra!

%% Organize data
if isempty(data(size(data,1),1))
    numsamples=size(data,1)-1; % avoid last empty line
else
    numsamples=size(data,1);
end
Samplename=data(2:numsamples,1);
signals=str2double(data(2:numsamples,6));
errors=str2double(data(2:numsamples,7));
if sum(isnan(errors))>0.5*numel(errors) % if no STD
    warning('Most of the Raw.STD empty! Using Raw.Average*Raw.RSD/100 instead...')
    errors=abs(signals.*str2double(data(2:numsamples,9)))/100;
end
nominal=str2double(data(2:numsamples,8));
analytename=data(2:numsamples,3);
elementname=data(2:numsamples,4);
Loc=str2double(data(2:numsamples,2));
numsamples=numsamples-1;

%% Select blanks
blanknames=[{'blk'},{'blank'}];
for n=1:numsamples
    if isnan(nominal(n))
        for m=1:length(blanknames)
            if ~isempty(regexpi(Samplename{n},blanknames{m})) % check if this line corresponds to a blank
                %                 disp([Samplename{n} ' is a blank (guess)'])
                nominal(n)=0;
            end
        end
    elseif nominal(n)==0
        %         disp([Samplename{n} ' is a blank (stated)'])
    end
end


%% Copy repeated STDs
for n=1:numsamples
    if isnan(nominal(n))
        sel=find(strcmp(Samplename{n}, Samplename) & ~isnan(nominal) & strcmp(analytename{n}, analytename),...
            1,'first'); % find a previous line with the same sample-name and analyte AND a nominal value
        if ~isempty(sel)
            nominal(n)=nominal(sel);
%             disp([Samplename{n} ' nominal data copied from ' Samplename{sel}])
        end
    end
end

%% Show and (de)select blnks and standards
selectednominals=~isnan(nominal);
list_strings=Samplename;
for n=1:numsamples
    if selectednominals(n)
    list_strings{n}=[list_strings{n} '-' analytename{n} '-' num2str(nominal(n)) 'ppm'];
    else
        list_strings{n}=[list_strings{n} '-' analytename{n} '-UNKNOWN'];
    end
end

% allow the user to select the nominal data to be used from now on
[selected,~] = listdlg('PromptString','Select blanks and datandards:',...
    'ListSize',[600 300],...
    'SelectionMode','multiple',...
    'ListString',list_strings,...
    'InitialValue',find(selectednominals));

% remove unwanted blanks and standards
for n=1:numsamples
    if sum(n==selected)==1
        if isnan(nominal(n)) % if lines with no data are selected, we will consider them as blanks!
            warning([Samplename{n} '-' analytename{n} ': no nominal data. It will be considered as a blank. If needed, consider including its nominal value in the  ' xstring ' column.'])
            nominal(n)=0;
        end
    else
        nominal(n)=NaN; % remove nominal 
    end
end

%% Calculate and plot calibration
STD=~isnan(nominal); % define standards and blanks
numanalytes=length(unique(analytename));
numelements=length(unique(elementname));
numuniquesamples=length(unique(Samplename));
rowsinplot=ceil(numanalytes^0.5);
rowsinplot2=ceil(numelements^0.5);

% create arrays
corrsignals=0.*signals;
correrrors=0.*signals;
totalerrors=0.*signals;
exterrors=0.*signals;

% start figure
figcalib = figure('units','normalized','outerposition',[0 0 1 1]);
hold on
step=0;

resultsstring=['Analyte' '	' 'Drift(total percent)' '	' 'formula' '	' 'R2' '	' 'Chisq'];

for analyte=unique(analytename)'
    step=step+1;
    subplot(rowsinplot,ceil(numanalytes/rowsinplot),step) % select plot position in the figure
    hold on % do not remove lines qhen we plot something new
    sel=(STD & strcmp(analyte, analytename) & ~isnan(signals)); % select standards
    sel2=(strcmp(analyte, analytename)); % select all measuremets for this analyte
    
    % standards only
    y=nominal(sel);
    x=signals(sel);
    t=Loc(sel); % location in the run (for drift calculations)
    dx=errors(sel);
    dy=0.01.*y; % consider 1% error in standards
    
    % all measurements
    x2=signals(sel2);
    dx2=errors(sel2);
    t2=Loc(sel2);
    
    
    if sum(y>0)>0 % if we have enough standards (not blanks)
        
        % calculate drift
        drift=linear_regression_chisq_fn(t(y>0),0.*t(y>0)+1,x(y>0)./y(y>0),((dx(y>0)./y(y>0)).^2+(dy(y>0).*x(y>0)./y(y>0).^2).^2).^0.5);
        totaldrifttemp=(drift.intercepts(1)+drift.slopes(1)*max(t2))/drift.intercepts(1)-1;
        if drift.sloping && abs(totaldrifttemp)<1 % if the linera regression is NOT compatible with a horizontal slope, and drift is less than 100% (>100% drift means that samething is VERY wrong)
            totaldrift=(drift.intercepts(1)+drift.slopes(1)*max(t2))/drift.intercepts(1)-1;
        else % if not
            totaldrift=0; % no drift
        end
        
        % drift correction
        x=x./(1+totaldrift*t/max(t2));
        x2=x2./(1+totaldrift*t2/max(t2));
        
        % calibration
        myfit=linear_regression_chisq_fn(x,dx,y,dy); % calculate the solopes and intercepts 
        
        y2=myfit.intercepts(1)+myfit.slopes(1)*(x2); % this is the data calibrated
        y2max=max(myfit.intercepts+myfit.slopes.*(x2'))'; % this is the positive calibration error
        y2min=min(myfit.intercepts+myfit.slopes.*(x2'))'; % this is the negative calibration error
        
        
        corrsignals(sel2)=y2;
        correrrors(sel2)=myfit.slopes(1)*dx2; % transmission of the measurement uncertainties (std). Internal error.
        totalerrors(sel2)=(correrrors(sel2).^2+((y2max-y2min)/2).^2).^0.5; % internal error + mean calibration error
        
        

        
        eq=[num2str(myfit.intercepts(1)) '+x*' num2str(myfit.slopes(1))]; % equation in text format
        eqerror=[' '];
        
        BEC(step)=(max(myfit.intercepts)-min(myfit.intercepts))/2; % blank equivalent concentration
        
        resultsstring=[resultsstring '	\n' analyte{:} '	' num2str(round(totaldrift*1000)/10) '	' 'y=' eq '	'];
        resultsstring=[resultsstring '	' num2str(myfit.R2) '	' num2str(myfit.redchisq) '	'];
        
        
         plot(y2,x2,'x','Color',[0.6 0.6 0.6]) % plot calibrated data
        
        % calculate calibration lines
        xplot=linspace(min(signals(sel2)-errors(sel2)),max(signals(sel2)+errors(sel2))*1.2,300);
        yplot=myfit.intercepts(1)+myfit.slopes(1)*(xplot);
        yplotmax=yplot;
        yplotmin=yplot;
        for n=1:length(myfit.intercepts)
            yploti=myfit.intercepts(n)+myfit.slopes(n)*(xplot);
            yplotmax=max(yplotmax,yploti);
            yplotmin=min(yplotmin,yploti);
        end
        
        % plot calibration and uncertainties
        plot(yplot,xplot,'-r')
        plot(yplotmax,xplot,'--r')
        plot(yplotmin,xplot,'--r')
        
        % plot standards
        for z=1:length(x)
            plot(y(z),x(z),'xb')
            plot([y(z)-dy(z),y(z)+dy(z)],[x(z),x(z)],'-b')
        end
        
        % modify plot limits (crop the graphs)
        xlim([0.00 max(corrsignals(sel2))*1.2])
        ylim([min(signals(sel2)) max(signals(sel2))*1.2])
        % 		xlim([min(abs(corrsignals(sel2)+correrrors(sel2))) max(corrsignals(sel2)+correrrors(sel2))*1.2])
        %         set(gca,'xscale','log')
        %
        % 		ylim([min(abs(signals(sel2)+errors(sel2))) max(signals(sel2)+errors(sel2))*1.2])
        %         set(gca,'yscale','log')
        
        % labels and titles
        xlabel('nominal')
           if step==1
               ylabel('ICP (drift corr.)')
           else
               ylabel('ICP')
           end
        
        grid on
        title(analyte)
        
    end
end

% %% export first calibration figure
% outputfilename=[profilefile '_calib_all'];
% fnam=strcat(outputfilename,'.eps');
% % fnam=fnam{1};
% snam='MSWord'; % note: NO extension...
% s=hgexport('readstyle',snam);
% hgexport(figcalib,fnam,s);
% % close(figcalib)

%% export next data to text file
clc
outputfilename=[profilefile '_resutls' '.txt'];
delete(outputfilename)
diary(outputfilename) % start writing the results file

disp(['**************************************'])
disp(['ICPOES data reduction'])
disp(['Version: ',mfilename])
disp(['Author: Angel Rodes (SUERC)'])
disp(['ICP data file: ' profilefile])
disp(['Date: ',datestr(clock,0)])
disp(['Paste this data in excel/calc'])
disp(['**************************************'])

%% print calibration
disp([' '])
disp(['---- Calibration ----'])
%     disp('Linear regressions')
%
disp(sprintf(resultsstring))

%% check standards

disp([' '])
disp(['Analyte' '	' 'BEC' '	' 'LOD' '	' 'LOQ' '	' 'stdev'])
count=0;
bec=BEC;
lod=3*BEC;
loq=10*BEC;
for ana=unique(analytename)'
    count=count+1;
    sel=(STD & nominal>0 & strcmp(ana, analytename) & ~isnan(corrsignals));
    scatterstd(count)=sum(abs(corrsignals(sel)-nominal(sel))./corrsignals(sel))/sum(sel);
    disp([ana{:} '	' num2str(bec(count),2) '	' num2str(lod(count),2) '	' num2str(loq(count),2) '	' num2str(scatterstd(count)*100,2) '%'])
end

% LOD by element

	disp([' '])
    disp(['Element' '	' 'LOD' '	' 'LOQ'])
	count=0;
	for ana=unique(elementname)'
		count=count+1;
		sel=(STD & nominal>0 & strcmp(ana, elementname));
		sel2=(STD & nominal==0 & ~isnan(corrsignals) & ~isnan(correrrors) & strcmp(ana, elementname));
		if sum(sel2)==1
			anamean=(corrsignals(sel2));
			anaerror2=(totalerrors(sel2));
            bec2(count)=max(anamean,0)+anaerror2;
			lod2(count)=max(anamean,0)+3*anaerror2;
			loq2(count)=max(anamean,0)+10*anaerror2;
		elseif sum(sel2)==0
			sel3=(STD & ~isnan(corrsignals) & ~isnan(correrrors) & strcmp(ana, elementname));
            bec2(count)=median(totalerrors(sel3));
			lod2(count)=median(totalerrors(sel3)*3);
			loq2(count)=median(totalerrors(sel3)*10);
		else
			anamean=(sum(corrsignals(sel2)./totalerrors(sel2).^2)/sum(1./totalerrors(sel2).^2));
			anaerror2=max(std(corrsignals(sel2)),median(totalerrors(sel2)));
            bec2(count)=max(anamean,0)+anaerror2;
			lod2(count)=max(anamean,0)+3*anaerror2;
			loq2(count)=max(anamean,0)+10*anaerror2;
		end
		disp([ana{:} '	' num2str(lod2(count),2) '	' num2str(loq2(count),2)])
	end

%% print results

disp(['-------------'])
disp(['Results line by line'])
disp(['-------------'])
disp(['Samplename' '	' 'analytename'  '	' 'signal'  '	' 'error'  '	' 'nominal'  '	' 'concentration'  '	' 'intrerror'  '	' 'totalrerror'])

for i=1:numsamples
    disp([Samplename{i} '	' analytename{i}  '	' num2str(signals(i))  '	' num2str(errors(i))  '	' num2str(nominal(i))  '	' num2str(corrsignals(i))  '	' num2str(correrrors(i)) '	' num2str(totalerrors(i))])
end



%% print table with all data


disp(['-------------'])
disp(['Tabulated by analyte'])
disp(['-------------'])

headstr=['Samplename' '	'];
for ana=unique(analytename)'
    headstr=[headstr ana{:} '	' ana{:} ' int. error	' ana{:} ' total error	'];
end
disp(headstr)

for location=unique(Loc)'
    sel=(Loc==location);
    selname=find(sel,1,'first');
    samplestr=[Samplename{selname} '	'];
    count=0;
    for ana=unique(analytename)'
        count=count+1;
        sel2=(sel & strcmp(ana, analytename) & ~isnan(corrsignals) & ~isnan(correrrors));
        if sum(sel2)==1
            anamean=num2str(corrsignals(sel2));
            anaerror1=num2str(correrrors(sel2));
            anaerror2=num2str(totalerrors(sel2));
        elseif sum(sel2)==0
            anamean=' ';
            anaerror1=' ';
            anaerror2=' ';
        else
            anamean=num2str(sum(corrsignals(sel2)./correrrors(sel2).^2)/sum(1./correrrors(sel2).^2));
            anaerror1=num2str(1/sum(1./correrrors(sel2).^2)^0.5);
            anaerror2=num2str(1/sum(1./totalerrors(sel2).^2)^0.5);
        end
        if anamean<lod(count)
            %             samplestr=[samplestr '-' '	' '-' '	' '-' '	']; % do not show <LOD
            samplestr=[samplestr anamean '	' anaerror1 '	' anaerror2 '	'];
        else
            samplestr=[samplestr anamean '	' anaerror1 '	' anaerror2 '	'];
        end
    end
    disp(samplestr)
    
    if sum(sel & STD)>0
        selname=find(sel,1,'first');
        STDstr=[Samplename{selname} ' NOMINAL	'];
        for ana=unique(analytename)'
            sel2=(sel & strcmp(ana, analytename) & ~isnan(nominal));
            if sum(sel2)==1
                anamean=num2str(nominal(sel2));
                anaerror1=' ';
                anaerror2=' ';
            elseif sum(sel2)==0
                anamean=' ';
                anaerror1=' ';
                anaerror2=' ';
            else
                anamean=num2str(mean(nominal(sel2)));
                anaerror1=num2str(std(nominal(sel2),1));
                anaerror2=num2str(std(nominal(sel2),1));
                if std(nominal(sel2),1)<mean(nominal(sel2))*0.00001+0.000001
                    anaerror1=' ';
                    anaerror2=' ';
                end
            end
            STDstr=[STDstr anamean '	' anaerror1 '	' anaerror2 '	'];
        end
        disp(STDstr)
    end
end



%% print "commercial" table tabulated
disp(['--------------------'])
disp(['Tabulated by element'])
disp(['--------------------'])

n=0;

headstr=['Samplename' '	'];
for ana=unique(elementname)'
    headstr=[headstr ana{:} '	±	' 'error	'];
end
disp(headstr)

countlocation=0;
for location=unique(Loc)' % for each location
    countlocation=countlocation+1;
    sel=(Loc==location);
    selname=find(sel,1,'first');
    %         if sum(sel & STD)==0
    samplestr=[Samplename{selname} '	'];
    count=0;
    for ana=unique(elementname)' % for each element
        count=count+1;
        sel2=(sel & strcmp(ana, elementname) & ~isnan(corrsignals) & ~isnan(correrrors));
        if sum(sel2)==1
            anamean=corrsignals(sel2);
            anaerror1=correrrors(sel2);
            anaerror2=totalerrors(sel2);
        elseif sum(sel2)==0
            anamean=NaN;
            anaerror1=NaN;
            anaerror2=NaN;
        else
            if min(correrrors(sel2))>0 % if we have uncertainties
                anamean=sum(corrsignals(sel2)./totalerrors(sel2).^2)/sum(1./totalerrors(sel2).^2); % weighted mean
                anaerror1=1/sum(1./correrrors(sel2).^2)^0.5; % standard error of the weighted mean (internal)
                anaerror2=1/sum(1./totalerrors(sel2).^2)^0.5; % standard error of the weighted mean (external)
            else % if we DO NOT have uncertainties for some reason
                anamean=mean(corrsignals(sel2)); % mean
                anaerror1=std(corrsignals(sel2),1); % standard deviation
                anaerror2=std(corrsignals(sel2),1);
            end
        end
        anaerror2=max(bec2(count),min(anaerror2,anamean)); % external uncertainties should be at least the blank equivalent concentration
        
        % get nominal value (if any)
        if sum(sel2)==0
            nomvalue=NaN;
        else
            nomvalue=nominal(find(sel2,1,'first'));
        end
        
        n=n+1;
        simpletable.location(n)=location;
        simpletable.element(n)=ana;
        simpletable.nominal(n)=nomvalue;
        simpletable.mean(n)=anamean;
        simpletable.error(n)=anaerror2;
        
        % save string
        precis=floor(log10(bec2(count))-1); % calculate precison to show only relevant figures based on BEC
        if sum(sel2)>0 % if we have data
            if anamean<lod2(count) % if below limit of detection
                %             samplestr=[samplestr '-' '	' '-' '	' '-' '	']; % do not
                %             show <LOD
                anastr=[ ' ' '	<	' num2str(round(lod2(count)/10^precis)*10^precis)];
            else % if above limit of detection
                anastr=[num2str(round(anamean/10^precis)*10^precis) '	±	' num2str(round(anaerror2/10^precis)*10^precis)];
            end
        else % if no data
            anastr=[ ' ' '	 	' ' '];
        end
        samplestr=[samplestr anastr '	']; % add data to string
        elemmean(count,countlocation)=anamean; % store mean
        elemerror(count,countlocation)=anaerror2; % store external error
        
    end
    disp(samplestr) % diplay string
    %         end
end


%% check elements
disp(['-------------'])
disp(['Check element averages'])
disp(['-------------'])
disp(['Element' '	' 'N' '	' 'SDOM' '	' 'mean error'])

figelements = figure('units','normalized','outerposition',[0 0 1 1]);
hold on
count=0;
xplot=linspace(0,2,1000);
% xplot=logspace(-2,2,1000);
for analyte=unique(elementname)'
    count=count+1;
    subplot(rowsinplot2,ceil(numelements/rowsinplot2),count)
    hold on
    yplot=0.*xplot;
    countlocation=0;
    checkcount=0;
    for location=unique(Loc)'
        countlocation=countlocation+1;
        sel=(Loc==location);
        sel2=(sel & strcmp(analyte, elementname) & ~isnan(corrsignals) & ~isnan(correrrors));
        mus=corrsignals(sel2)/elemmean(count,countlocation);
        sigmas=totalerrors(sel2)/elemmean(count,countlocation);
        for j=1:length(mus)
            if ~isnan(sigmas(j))
                if elemmean(count,countlocation)>lod2(count)
                    yplot=yplot+normpdf_local(xplot,mus(j),sigmas(j));
                    if isnan(max(normpdf_local(xplot,mus(j),sigmas(j))))
                        disp([num2str(mus(j)) ' +/- ' num2str(sigmas(j))])
                    end
                    checkcount=checkcount+1;
                end
            end
        end
    end
    if max(yplot)>0
        area(xplot,yplot/sum(yplot))
        
        % genereate numbers based on (xplot,yplot) pdf
        nrandnumbers=10000;
        nout=round(yplot/sum(yplot)*nrandnumbers);
        accumnout=cumsum(nout);
        for k=1:sum(nout)-1
            pos(k)=find(accumnout>k,1,'first');
        end
        randomnumbersx=xplot(pos);
        
        
        meanerrorelememt=1;
        refcurve=yplot/sum(yplot);
        diffcurves=sum(refcurve);
        for err=0.0005:0.0005:1
            curve1=normpdf_local(xplot,1,err)/sum(normpdf_local(xplot,1,err));
            difference=sum(abs(refcurve-curve1));
            if difference<diffcurves
                meanerrorelememt=err;
                diffcurves=difference;
            end
        end
        
        sdomelement=std(randomnumbersx)/checkcount^0.5;
        
        plot(xplot,normpdf_local(xplot,1,meanerrorelememt)/sum(normpdf_local(xplot,1,meanerrorelememt)),'--r')
        
        xlim([min(xplot(yplot>max(yplot)/100)) max(xplot(yplot>max(yplot)/100))])
        xlabel(['meas/mean = ' num2str(round(mean(randomnumbersx)*100)) ' ± ' num2str(round(meanerrorelememt*100)) ' %'])
        ylabel('P')
        title([analyte{:} ' (n=' num2str(checkcount) ')'])
        
        if max(yplot)==0
            text(1,0.5,'no data')
        end
        if isnan(yplot)
            text(0.8,0.5,'NaN data')
        end
        
        
        disp([analyte{:} '	' num2str(checkcount) '	' num2str(round(sdomelement*100*10)/10) '%' '	' num2str(round(meanerrorelememt*100*10)/10) '%'])
    end
end

% %% export figelements
% outputfilename=[profilefile '_figelements'];
% fnam=strcat(outputfilename,'.eps');
% % fnam=fnam{1};
% snam='MSWord'; % note: NO extension...
% s=hgexport('readstyle',snam);
% hgexport(figelements,fnam,s);
% % close(figcalib)

%% check drift
figdrift = figure('units','normalized','outerposition',[0 0 1 1]);
hold on
step=0;
for analyte=unique(analytename)'
	step=step+1;
	subplot(rowsinplot,ceil(numanalytes/rowsinplot),step)
	hold on
	sel=(STD & strcmp(analyte, analytename) & nominal>0  & ~isnan(corrsignals) & ~isnan(correrrors));

%     drift=fit(find(sel),corrsignals(sel)./nominal(sel),'poly1');
% %     plot(find(sel),drift(find(sel))./mean(drift(find(sel))),'--r')
%     plot(find(sel),find(sel).*0+1,'-g')
% %     plot(find(sel),corrsignals(sel)./nominal(sel),'xb')

    for j=find(sel)'
        scol=nominal(j)/max(nominal(sel));
        symbolcolor=[min(1,scol*2) 0 min(1,(1-scol)*2)]; % from blue to red
        plot(j,corrsignals(j)./nominal(j),'xb','Color',symbolcolor)
        plot([j,j],[(corrsignals(j)+totalerrors(j))/nominal(j),(corrsignals(j)-totalerrors(j))/nominal(j)],'-b','Color',symbolcolor)
    end

    meanerror=mean(totalerrors(sel)./corrsignals(sel));

%     ylim([1-meanerror*2 1+meanerror*2])
    xlabel('analysis')
	ylabel('concentrations/nominal')
	grid on
	title(analyte)
    xlim([0 length(STD)])

end
%
% %% export figdrift
% outputfilename=[profilefile '_drift'];
% fnam=strcat(outputfilename,'.eps');
% % fnam=fnam{1};
% snam='MSWord'; % note: NO extension...
% s=hgexport('readstyle',snam);
% hgexport(figdrift,fnam,s);
% % close(figcalib)

%% check drift2
figdrift = figure('units','normalized','outerposition',[0 0 1 1]);
hold on
step=0;
for analyte=unique(analytename)'
    step=step+1;
    subplot(rowsinplot,ceil(numanalytes/rowsinplot),step)
    hold on
    sel=(STD & strcmp(analyte, analytename)  & ~isnan(corrsignals) & ~isnan(correrrors));
    
    plot(find(sel),find(sel).*0,'-g')
    plot(find(sel),corrsignals(sel)-nominal(sel),':r')
    
    for j=find(sel)'
        scol=nominal(j)/max(nominal(sel));
        symbolcolor=[min(1,scol*2) 0 min(1,(1-scol)*2)]; % from blue to red
        plot(j,corrsignals(j)-nominal(j),'xb','Color',symbolcolor)
        plot([j,j],[(corrsignals(j)+totalerrors(j))-nominal(j),(corrsignals(j)-totalerrors(j))-nominal(j)],'-b','Color',symbolcolor)
    end
    
    
    xlabel('analysis')
    if step==1
        ylabel('concentration-nominal')
    else
        ylabel('dev. ppm')
    end
    grid on
    title(analyte)
    xlim([0 length(STD)])
    sel3=(STD & ~isnan(corrsignals) & ~isnan(correrrors));
    maxerrorplot=prctile(abs(corrsignals(sel3)-nominal(sel3)),99);
    ylim([-maxerrorplot +maxerrorplot])
    
end

% %% export figdrift2
% outputfilename=[profilefile '_drift2'];
% fnam=strcat(outputfilename,'.eps');
% % fnam=fnam{1};
% snam='MSWord'; % note: NO extension...
% s=hgexport('readstyle',snam);
% hgexport(figdrift,fnam,s);
% % close(figcalib)

%% check errors
figdrift = figure('units','normalized','outerposition',[0 0 1 1]);
hold on
step=0;
for analyte=unique(analytename)'
    step=step+1;
    subplot(rowsinplot,ceil(numanalytes/rowsinplot),step)
    hold on
    sel=(STD & strcmp(analyte, analytename)  & ~isnan(corrsignals) & ~isnan(correrrors));
    
    
    xdata=corrsignals(sel);
    ydata=abs(corrsignals(sel)-nominal(sel));
    for n=1:20
        xdata=[xdata;corrsignals(sel)];
        ydata=[ydata;abs(normrnd_BoxMuller(corrsignals(sel),totalerrors((sel)))-nominal(sel))];
    end
    ydata=ydata;
%     drift=fit(xdata,ydata,'poly1');
%     plot(sort(xdata),drift(sort(xdata)),'--g')
%     plot(sort(xdata),-drift(sort(xdata)),'--g')
    
    plot(sort(xdata),sort(xdata).*0,'-g')
    
    for j=find(sel)'
        scol=nominal(j)/max(nominal(sel));
        symbolcolor=[min(1,scol*2) 0 min(1,(1-scol)*2)]; % from blue to red
        plot(corrsignals(j),(corrsignals(j)-nominal(j)),'xb','Color',symbolcolor)
        plot([corrsignals(j),corrsignals(j)],[((corrsignals(j)+totalerrors(j))-nominal(j)),((corrsignals(j)-totalerrors(j)))-nominal(j)],'-b','Color',symbolcolor)
    end
    
    
    xlabel('conc. ppm')
    if step==1
        ylabel('concentration-nominal')
    else
        ylabel('dev. ppm')
    end
    grid on
    title(analyte)
    %     xlim([0 length(STD)])
    sel3=(STD & ~isnan(corrsignals) & ~isnan(correrrors));
    maxerrorplot=prctile(abs(corrsignals(sel3)-nominal(sel3))+totalerrors(sel3),95);
    ylim([-maxerrorplot +maxerrorplot])
    
end

% %% export errors
% outputfilename=[profilefile '_errors_vs_conc'];
% fnam=strcat(outputfilename,'.eps');
% % fnam=fnam{1};
% snam='MSWord'; % note: NO extension...
% s=hgexport('readstyle',snam);
% hgexport(figdrift,fnam,s);
% % close(figcalib)


%% bye!
disp('---end of text output---')
diary off % stop writting results file