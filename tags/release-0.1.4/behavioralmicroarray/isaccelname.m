function v = isaccelname(prop)

accelnames = {'accmag','d2theta','absd2theta','abssmoothd2theta'};
v = ismember(prop,accelnames);