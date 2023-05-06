function handles = grey2binarySliderUpdater(handles,currentSliderValue)

    adjustedSliderValue  = currentSliderValue ;

	handles.constants.grey2BWthreshold = round(1e3*adjustedSliderValue)/1e3 ;

	set(handles.grey2binaryThresholdSlider,'Value', round(1e3*currentSliderValue)/1e3 )
	set(handles.grey2binaryThresholdSliderValue,'String',num2str( round(1e3*adjustedSliderValue)/1e3 ))

	if handles.imageLoaded
		set(handles.segmentationButton,'Enable','on') ;
	end
	
end

