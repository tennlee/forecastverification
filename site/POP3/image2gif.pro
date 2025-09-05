pro image2gif,giffile

;Reads image from the window and write to designated gif file
;This seems to be necessary when running IDL on the PC with true color

;Get color scale
tvlct,r,g,b,/get

;Read image in true color and convert to pseudocolor
image=color_quan(tvrd(true=1),1,r,g,b)

;Write to file
write_gif,giffile,image,r,g,b

end