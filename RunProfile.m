function profileInfo = RunProfile()
% RunProfile - Profile the full bubble-detection pipeline in main.m.
%
% Usage:
%   profileInfo = RunProfile();
%
% Results are saved in ProfileResults/profileInfo.mat and an HTML profiler
% report is generated in ProfileResults/html.

outputFolder = 'ProfileResults';
htmlFolder = fullfile(outputFolder, 'html');

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

if ~exist(htmlFolder, 'dir')
    mkdir(htmlFolder);
end

profile off;
profile clear;
profile on;

try
    runMainScript();
catch caughtError
    profile off;
    profileInfo = profile('info');
    save(fullfile(outputFolder, 'profileInfo_failed.mat'), 'profileInfo');
    rethrow(caughtError);
end

profile off;
profileInfo = profile('info');
save(fullfile(outputFolder, 'profileInfo.mat'), 'profileInfo');
profsave(profileInfo, htmlFolder);
profile viewer;

fprintf('Profile saved to %s\n', fullfile(outputFolder, 'profileInfo.mat'));
fprintf('HTML profile report saved to %s\n', htmlFolder);

end

function runMainScript()
run('main.m');
end
