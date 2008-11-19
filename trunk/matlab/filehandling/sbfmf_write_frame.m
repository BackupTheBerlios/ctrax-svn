function sbfmf_write_frame(fp,stamp,idx,val)

npixels = length(idx);
fwrite(fp,npixels,'uint32');
fwrite(fp,stamp,'double');
fwrite(fp,uint32(idx-1),'uint32');
fwrite(fp,uint8(val),'uint8');
