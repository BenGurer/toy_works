function range = findRange(data)

ampMin = realmax;
ampMax = 0;
nScans = length(data);
for scan=1:nScans
  if ~isempty(data{scan})
     thisData = data{scan}(:);
     thisData = thisData(~isinf(thisData));
    ampMin = min([ampMin min(thisData)]);
    ampMax = max([ampMax max(thisData)]);
  end
end
if (ampMin <= ampMax)
  range = [ampMin ampMax];
else
  % if amp data is empty, need to make sure min < max
  range = [0 1];
end