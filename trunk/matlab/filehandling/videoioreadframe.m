function im = videoioreadframe(readerobj,f)

seek(readerobj,f);
im = getframe(readerobj);