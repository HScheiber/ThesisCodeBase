function nstcouple = Get_nstcouple(tau,dt)
% Function to generate an appropriate nstpcouple or nsttcouple value for a
% given timestep + tau_p or timestep + tau_t combo.
% nstcouple is the number of steps between coupling the thermostat or
% barostat, it should be an integer >= 1.
% Tau is the time constant, dt is the timestep
nstcouple = max(round(tau/(20*dt)),1);

end