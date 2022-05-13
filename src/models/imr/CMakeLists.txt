add_library(fabm_models_imr OBJECT
	imr_model_library.F90
	norwecom/norwecom.F90
	)

add_dependencies(fabm_models_imr fabm_base)
