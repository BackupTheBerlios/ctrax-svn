function trx = prune_short_trajectories(trx,minnframes)

keepidx = [trx.nframes] >= minnframes;
if ~all(keepidx),
  fprintf(['Ignoring flies [',num2str(find(~keepidx)),'], trajectory lengths < %d\n'],minnframes);
  trx = trx(keepidx);
end