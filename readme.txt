Matlab, choose "HOME"->"Set Path"
delete original "hdrvdp" and subpath from the Matlab path
add new "hdrvdp" and subpath in workspace to the Matlab path again(right click the "hdrvdp",choose "Add to Path")

Update 20170928
1, add JND method.
2, add only ssim_map, current the map is calculated by cropping rectangle.
3, add time feedback, min, max value of ssim map representation which could help setting the min and max hdrvdp parameters.
4, add save button to save the result image.

Update 20171009
1, add vdm map
2, fix save function for cortexMask result
3, fix bottom-up selection and top-down selection 

Update 20171010
1, using the same color map to hdrvdp and vdm
2, add color scale image
3, move the text result box
4, add “PPD” parameter, means pixels per degree

Update 20171019
1, add jnd value for ssim/cwssim/vdm in result
2, change the ssim/cwssim index map
3, fix the error when WNum get 40

Update 20180113
1, add SaltencyMap function and insert it to tools.
2, replace SaltencyMap to SSim