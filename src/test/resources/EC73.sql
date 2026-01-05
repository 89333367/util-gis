select
    did,
    jobStartTime,
    jobEndTime,
    jobArea,
    effectiveJobArea,
    jobWidth
from
    farm_work
where
    jobArea>0.5
        and checkFlag = 1
        and year(jobEndTime)= 2025
  and MONTH(jobEndTime)= 11
  and did in('EC71BD2408280034'
    , 'EC71BD2408280043'
    , 'EC73BD2503190283'
    , 'EC73BD2508220088'
    , 'EC73BD2504110747'
    , 'EC73BD2504110136'
    , 'EC73BD2508220108'
    , 'EC73BD2508220106'
    , 'EC73BD2509060269'
    , 'EC73BD2509060268'
    , 'EC73BD2508220107'
    , 'EC73BD2504110757'
    , 'EC73BD2504270161')
order by
    jobArea desc
