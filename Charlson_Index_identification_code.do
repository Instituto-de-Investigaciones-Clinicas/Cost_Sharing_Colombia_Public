*********************************************************************************************************
* 		Important file paths
*********************************************************************************************************
	*Log_files
	global logs_cs /home/log_files................
	*File path data bases CS:
	global bases_cs /home/data_bases..................
	*File path results
	global resultados_cs /home/results............
	
************************************************************************************************************
** Cost-Sharing in Medical Care Can Increase Adult Mortality: Evidence from Colombia
** Aim: Charlson index identification using monthly files of the UPC database from the ministry of health 

//  task:      	use the monthly UPC databases (administrative database) to identiy the disease included in the Charlson index
//  authors:    Giancarlo Buitrago, Javier Amaya-Nieto, Marcos Vera-Hernandez and Grant Miller 
//  date: 		February 20 / 2024
************************************************************************************************************

** Opening the log file
	log using "${logs_cs}/Charlson index identification.smcl", append


		foreach year in 2011(1)2019 {
		foreach month in 01 02 03 04 05 06 07 08 09 10 11 12 {
		use PersonaBasicaID DiagnosticoCD FechaServicio ambitosprocedimientocd ProcedimientoCD diasestancia tipocodigoprocedimentocd ValorPagado using "${bases_ori}/UPC/UPC`year'/upc_`year'_`month'.dta", clear

	******************************************************************************************
	** Charlson index identification code
	**-+-\=========================================***------------------
		gen cie10_3d= DiagnosticoCD
		recast str3 cie10_3d, force
		gen cie10_2d=DiagnosticoCD
		recast str2 cie10_2d, force

		gen iam=1 if cie10_3d=="I21" | cie10_3d=="I22" | DiagnosticoCD=="I249" | DiagnosticoCD=="I252" | DiagnosticoCD=="I256"
		gen icc=1 if cie10_3d=="I50" | DiagnosticoCD=="I130" 
		gen evp=1 if cie10_3d=="I71" | DiagnosticoCD=="I790" | DiagnosticoCD=="I739" | cie10_3d=="R02" | DiagnosticoCD=="Z958" | DiagnosticoCD=="Z959"

		gen acv=1 if cie10_3d=="I60" | cie10_3d=="I61" | cie10_3d=="I62" | cie10_3d=="I63" | cie10_3d=="I65" | cie10_3d=="I66" | DiagnosticoCD=="G450" /*
		*/ | DiagnosticoCD=="G451" | DiagnosticoCD=="G452" | DiagnosticoCD=="G458" | DiagnosticoCD=="G459" | cie10_3d=="G46" | cie10_3d=="I64" /*
		*/ | DiagnosticoCD=="G454" | DiagnosticoCD=="I670" | DiagnosticoCD=="I671" | DiagnosticoCD=="I672" | DiagnosticoCD=="I674" | DiagnosticoCD=="I675" /*
		*/ | DiagnosticoCD=="I676" | DiagnosticoCD=="I677" | DiagnosticoCD=="I678" | DiagnosticoCD=="I679" | DiagnosticoCD=="I681" | DiagnosticoCD=="I682" /*
		*/ | DiagnosticoCD=="I688" | cie10_3d=="I69"

		gen deme=1 if cie10_3d=="F00" | cie10_3d=="F01" | cie10_3d=="F02" | DiagnosticoCD=="F051"

		gen epc=1 if cie10_3d=="J40" | cie10_3d=="J41" | cie10_3d=="J42" | cie10_3d=="J44" | cie10_3d=="J43" | cie10_3d=="J45" | cie10_3d=="J46" /*
		*/ | cie10_3d=="J47" | cie10_3d=="J67" | cie10_3d=="J44" | cie10_3d=="J60" | cie10_3d=="J61" | cie10_3d=="J62" | cie10_3d=="J63" | cie10_3d=="J66" /*
		*/ | cie10_3d=="J64" | cie10_3d=="J65"

		gen enf_conec=1 if cie10_3d=="M32" | cie10_3d=="M34" | DiagnosticoCD=="M332" | DiagnosticoCD=="M053" | DiagnosticoCD=="M058" | DiagnosticoCD=="M059" /*
		*/ | DiagnosticoCD=="M060" | DiagnosticoCD=="M063" | DiagnosticoCD=="M069" | DiagnosticoCD=="M050" | DiagnosticoCD=="M052" | DiagnosticoCD=="M051" /*
		*/ | DiagnosticoCD=="M353"

		gen ulc_pep=1 if cie10_3d=="K25" | cie10_3d=="K26" | cie10_3d=="K27" | cie10_3d=="K28"

		gen enf_hep=1 if DiagnosticoCD=="K702" | DiagnosticoCD=="K703" | cie10_3d=="K73" | DiagnosticoCD=="K717" | DiagnosticoCD=="K740" | DiagnosticoCD=="K742" /*
		*/ | DiagnosticoCD=="K746" | DiagnosticoCD=="K743" | DiagnosticoCD=="K744" | DiagnosticoCD=="K745"

		gen diab=1 if DiagnosticoCD=="E109" | DiagnosticoCD=="E119" | DiagnosticoCD=="E139" | DiagnosticoCD=="E149" | DiagnosticoCD=="E101" | DiagnosticoCD=="E111" /*
		*/ | DiagnosticoCD=="E131" | DiagnosticoCD=="E141" | DiagnosticoCD=="E105" | DiagnosticoCD=="E115" | DiagnosticoCD=="E135" | DiagnosticoCD=="E145"

		gen diab_compl=1 if DiagnosticoCD=="E102" | DiagnosticoCD=="E112" | DiagnosticoCD=="E132" | DiagnosticoCD=="E142" | DiagnosticoCD=="E113" /*
		*/ | DiagnosticoCD=="E133" | DiagnosticoCD=="E143 E104" | DiagnosticoCD=="E114" | DiagnosticoCD=="E134" | DiagnosticoCD=="E144" | DiagnosticoCD=="E103"

		gen paraple=1 if cie10_3d=="G81" | DiagnosticoCD=="G041" | DiagnosticoCD=="G820" | DiagnosticoCD=="G821" | DiagnosticoCD=="G822"

		gen enf_ren=1 if cie10_3d=="N03" | DiagnosticoCD=="N052" | DiagnosticoCD=="N053" | DiagnosticoCD=="N054" | DiagnosticoCD=="N055" | DiagnosticoCD=="N056" /*
		*/ | DiagnosticoCD=="N072" | DiagnosticoCD=="N073" | DiagnosticoCD=="N074" | cie10_3d=="N01" | cie10_3d=="N18" | cie10_3d=="N19" | cie10_3d=="N25" | DiagnosticoCD=="I120" /*
		*/ | DiagnosticoCD=="I131" | DiagnosticoCD=="N083"

		gen cancer=1 if cie10_2d=="C0" | cie10_2d=="C1" | cie10_2d=="C2" | cie10_2d=="C3" | cie10_3d=="C40" | cie10_3d=="C41" | cie10_3d=="C43" | cie10_3d=="C45" /*
		*/ | cie10_3d=="C46" | cie10_3d=="C47" | cie10_3d=="C48" | cie10_3d=="C49" | cie10_2d=="C5" | cie10_2d=="C6" | cie10_3d=="C70" | cie10_3d=="C71" /*
		*/ | cie10_3d=="C72" | cie10_3d=="C73" | cie10_3d=="C74" | cie10_3d=="C75" | cie10_3d=="C76" | cie10_3d=="C80" | cie10_3d=="C81" | cie10_3d=="C82" /*
		*/ | cie10_3d=="C83" | cie10_3d=="C84" | cie10_3d=="C85" | cie10_3d=="C883" | DiagnosticoCD=="C887" | DiagnosticoCD=="C889" | DiagnosticoCD=="C900" /*
		*/ | DiagnosticoCD=="C901" | cie10_3d=="C91" | cie10_3d=="C92" | cie10_3d=="C93" | cie10_3d=="C940" | DiagnosticoCD=="C941" | DiagnosticoCD=="C942" /*
		*/ | DiagnosticoCD=="C943" | DiagnosticoCD=="C9451" | DiagnosticoCD=="C947" | cie10_3d=="C95" | cie10_3d=="C96"

		gen can_meta=1 if cie10_3d=="C77" | cie10_3d=="C78" | cie10_3d=="C79" | cie10_3d=="C80"

		gen enf_hep_seve=1 if DiagnosticoCD=="K729" | DiagnosticoCD=="K766" | DiagnosticoCD=="K767" | DiagnosticoCD=="K721"

		gen vih=1 if cie10_3d=="B20" | cie10_3d=="B21" | cie10_3d=="B22" | cie10_3d=="B23" | cie10_3d=="B24"

		gen hta=1 if cie10_3d=="I10" | cie10_3d=="I11" | cie10_3d=="I12" | cie10_3d=="I13" | cie10_3d=="I14" | cie10_3d=="I15" 
		
	*******************************************************************************************************************
	** Colapsing data to obtain one register per person per month 
	
	//1. Collapse by person (the code runs for each month)
	collapse (max) iam icc evp acv deme epc enf_conec ulc_pep enf_hep diab diab_compl paraple enf_ren cancer can_meta enf_hep_seve vih hta, by(PersonaBasicaID)
	

	** replacing ceros 
		foreach var in iam icc evp acv deme epc enf_conec ulc_pep enf_hep diab diab_compl paraple enf_ren cancer can_meta enf_hep_seve vih hta {
			replace `var'=0 if `var'==.
			}
			
		** Charlson index generation 
		gen ind_char=hta+iam+icc+evp+acv+dem+epc+enf_conec+ulc_pep+enf_hep+diab+(diab_compl*2)+(paraple*2)+(enf_ren*2)+(cancer*2)+(can_meta*3)+(enf_hep_sev*3)+(vih*6)
	
		
	** Creating variables of year and month to identity where observations comes from
		gen month=`month'
		gen year=`year'
		
		
	** Saving databases
		compress 
		
		save "${bases_cs}/UPC/healt-serv_`year'_`month'.dta", replace
		}
		}
		
*****************************************************************************************************************************************
log close
