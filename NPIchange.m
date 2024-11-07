%% Function that work on data for NPI restrictions

%% ***************** INPUT ******************
% Mix_H - Mixing in the house.
% Mix_S - Mixing in school.
% Mix_W - Mixing iat workplaces.
% Mix_O - Mixing in other places.

%% ************** OUTPUT **********************

% ageMix - Age mixing metrice
function[ageMix] = NPIchange(Mix_H,Mix_S,Mix_W,Mix_O,npi)

if npi == "TOTAL_LOCKDOWN"
    ageMix = 1.5*Mix_H + Mix_W*0.3;
end

if npi == "LOCKDOWN_SCHOOL"
    ageMix = Mix_H + Mix_W*0.6;
end

if npi == "EASING"
    ageMix = 0.5*Mix_W + Mix_H + 0.3*Mix_O;
end

if npi == "RELAXED_RESTRICTION"
    ageMix = 0.5*Mix_W + Mix_H + 0.3*Mix_O+ Mix_S;
end

if npi == "RELAXED_RESTRICTION_SCHOOL"
    ageMix = 0.9*Mix_W + Mix_H + Mix_O + Mix_S;
end


end

