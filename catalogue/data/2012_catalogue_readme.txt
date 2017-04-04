2011-12-16

*******************************************************************************************************************

This text describes the development of the earthquake databse that will form the basis of the 2011 Australian Earthquake Hazard Map.  The base directory used for this development is: N:\ehp\georisk_earthquake\hazard_risk\Catalogues\

The catalogue is a composite catalogue and draws on information from the following existing catalogues:

1) GG-Cat:  An earthquake catalogue compiled by Gary Gibson.
   - Covers the area 110/156/-48/-10 from 1788-06-22 to 2010-05-26

2) QUAKES:  GA's catalogue of Australian and regional events.  
   - Covers the area 110/155/-45/-9 from 1902-05-07 to 2010-08-26

3) ISC AUST: All earthquakes in the International Seismological Centre's catalogue attributed to the network AUST.  This catalogue was primarily included to assist with the validation of catalogue magnitudes to determine the stations used to locate each event.
   - Covers the areas 111.856/155.235/-44.920/-10.350 from 1967-01-31 to 2008/04/30
   
4) ISC Timor/Banda Sea:  Additional earthquakes not captured in GG-Cat.  These data were added to calibrate recurrence for offshore seismic sources that may affect northern Australia.  Covers the area (from 1906/06/14 to 2011/04/17):
	- West: 108/110/-50/-10
	- North: 108/160/-10/-4
	- East: 156/160/-50/-10
	
5) QEDB: Queensland earthquake catalogue compiled by J. Rynn and D. Weatherly.  Covers state of Queensland and adjacent regions from 1866/12/29 to 2009/12/31
   
******************************************************************************************************************* 

Steps to create Master Catalogue (mdat):

1) Parse each of the individual catalogues into Matlab mat files using m-codes in respective directories.  For offshore earthquakes, use ISC data in "ISC_new" directory.  Run codes "parseISC" then "combineISC".  When parsing ISC locations with multiple location solutions, use the following logic to select preferred location:
		- if available, use EHB location 
		- select solution with lowest error with gap angle less than 180 degrees
		- if no solution applicable, select solution with maximum ISC solution ID assuming it is the most recent and robust
		- IDC preferred locations are removed for events of ISC median magnitude LT 5.0

2) Append GG-Cat to ISC (offshore) catalogue using "Merge ISC and GGCat\append_ISC_GG.m" and create the master catalogue "mdat"

3) Merge ISC AUST to mdat using "Merge MC and ISC AUST\merge_all_GG_ISC.m" for information on stations

4) Merge QUAKES to mdat using "Merge MC and QUAKES\merge_MC_QUAKES.m"

5) Merge QEDB to mdat using "Merge MC and QEDB\merge_MC_QEDB.m"

6) Create preferred catalogue parameters using "Preferred\getPrefMags.m"
		Logic for finding preferred source for each magnitude type:
		- MW: GG-Cat, HRVD, NEIC, AUST, Allen (2006, 2007), QEDB, Other
		- MS: GG-Cat, ISC, PAS, AUST, QEDB, Other
		- mb: GG-Cat, ISC, AUST, NEIC, IDC, QEDB, Other
		- ML: GG-Cat (MEL for lat <= -36 & lon >= 141 or lat <= -29 & lon >= 149 [from 1993 SCA network]), AUST, ISC, GG-Cat (MEL other), QEDB, Other
		- Other M: GG-Cat, QEDB
		
		* Preferred magnitude types as listed in order above (larger of mb/MS chosen)
		** Assume those earthquakes with unknown magnitude types are equivalent to ML

******************************************************************************************************************* 

Develop generic relations to correct magnitudes for unknown types that do not hav a RevML field prior to 1986

		Run "Preferred\get_M_corr_fact.m"
		
******************************************************************************************************************* 

Steps to revise magnitudes (see "Preferred\reviseMLs.m"):

		- Get approximate seismic station open and close times from first and last appearance in ISC catalogue
		- Assign zone for correcting magnitudes from the following polygons: 
			     walon = [129  110. 110.  135.0  135.0  138.3  138.3  129];
					 walat = [-10 -18.5  -45   -45   -29   -29   -10   -10];
					 ealon = [138.3  138.3  141.0  141.0  155.5  155.5 145.5 138.3];
					 ealat = [-10   -29   -29   -45   -45  -18 -10   -10];
			     salon = [135   135   141   141   135];                    
			     salat = [-29   -40   -40   -29   -29];  
			          
		- Get station list that recorded each earthquake, either from:
			1) ISC catalogue directly
			2) approximating stations within 1500 km of the epicentre assuming station open and close times from above
				a) If authority not MEL
						- If stations open from mid-2007, assume they currently open
						- If stations open prior to 1970, assume they were open from 1945 (estimated from GG-Cat)
						
				b) If authority MEL
						- Use SGData lookup table of MEL stations
			
		- With station distribution, undertake ML correction asumming methodology of Allen 2010 (AEES) under the following criteria
			
			1) Zone 1: Central & Western Australia
				a) Correct pre-1990 events to Gaull & Gregson 1990 assuming Richter (not ADE solutions)
				b) Correct pre-1968 ADE events to Gaull & Gregson 1990 assuming Richter
				c) Correct events with unknown magnitude types (inc MP, MD, MLI) to Gaull & Gregson 1990 assuming Richter
			
			2) Zone 2: Eastern Australia
				a) Correct pre-1990 events to Michael-Leiba & Manafant assuming Richter (not MEL solutions)
				b) Correct pre-1994 MEL events to Michael-Leiba & Manafant assuming Richter
				c) Correct post-2002 MEL events to Michael-Leiba & Manafant assuming Richter
				d) Correct 1994-2002 MEL events to Michael-Leiba & Manafant assuming Wilkie 96
				e) Correct events with unknown magnitude types (inc MP, MD, MLI) to Michael-Leiba & Manafant assuming Richter
				
			3) Zone 3: Mt Lofty & Flinders Ranges
				a) Correct pre-1968 ADE events to Greenhalgh & Singh 1986 assuming Richter
				b) Correct pre-1986 non-ADE events to Greenhalgh & Singh 1986 assuming Richter
				c) Correct events with unknown magnitude types (inc MP, MD, MLI) to Greenhalgh & Singh 1986 assuming Richter
				
			4) Apply generic magnitude correction for those earthquakes where no station information was availabale using coefficients determined using "get_M_corr_fact.m" which must be run after "reviseMLs.m".  The equation as of 2011-12-23 is: MLrev = 0.93*ML(Richter) + 0.02. Once this equation is updated, re-run "reviseMLs.m" (a bit circular, but necessary).

******************************************************************************************************************* 

Develop two ML to MW conversions for:
		1) central & western Australia (ML2MW\make_ML2MW_conv_WA.m)
		2) eastern Australia (ML2MW\make_ML2MW_conv_EA.m)
		
		Two seperate datasets were compiled of all known Mw estimates for Australian earthquakes, which also have ML estimates (WA_earthquakes_with_MW.txt & EA_earthquakes_with_MW.txt)
		
		Note: these codes add T. Allen Mw's to mdat_pref.MDAT_prefMW field where no existing entry exists.  Some ML's are also modified for key events based on professional opinion.  ML values are taken from GG-Cat.  The events are:
		
		ML 5.6 1989-12-27 23:24 (Newcastle - ML value from GG-Cat)
		ML 5.1 1994-08-06 11:03 (Ellalong - ML value from GG-Cat)
		ML 5.5 2011-04-16 5:31 (Bowen - ML value reviewed by T Allen)

*******************************************************************************************************************

Convert all magnitude types to MW and select preferred MW using "Preferred\get_pref_MW.m"

	1) Convert MS to MW using Scordilis (2006) for 3.0 <= MS <= 8.2 and h < 70 km
	2) Convert mb to MW using Scordilis (2006) for 3.5 <= mb <= 6.2
	3) Extrapolate Scordilis (2006) mb to MW conversion for AUST & MGO 3.0 <= mb < 3.5
	4) Convert ML to MW using Allen conversions for events within Australia (CWA, EA & MLFR)
	5) Now use Allen SEA - Convert ML to MW using Grunthal et al (2009) conversions for events outside Australia (i.e., additional ISC events)
	
	Select preferred MW based on following logic
	1) Conserve actual MW measurements
	2) Use MW from larger of MS/mb for Mx >= 6.0
	3) Use MW from revised ML (all magnitudes pre 1990)
	4) Use MW from ML (all magnitudes post 1990)
	5) Use MW from larger of MS/mb for Mx < 6.0
	
	Preferred MW catalogue, AUSTCAT.MW.VX.X.csv

	Select preferred M(existing) based on following logic with rescaling for ML
	1) Conserve actual MW measurements
	2) Use M from larger of MS/mb for Mx >= 6.0
	3) Use M from revised ML (all magnitudes pre 1990)
	4) Use M from ML (all magnitudes post 1990)
	5) Use M from larger of MS/mb for Mx < 6.0

******************************************************************************************************************* 
Versions:

V0.9:		Removes ISC mb's from outside Australia and uses Allen ML-MW conversions
V0.10: 	Replaces Grunthal ML conversion with Allen SEA conversion
V0.11:	- Adds a new column with Grunthal ML conversion (ML2MWG)
				- Replaces Scordilis (2006) mb2MW conversion with Allen Australian mb2MW conversion
				- If preferred MW comes from ML, uses ML2MWG
				- Revised date ranges for ML corrections based on discussions at EQ Hazard Workshop (2012-02-03)

******************************************************************************************************************* 
AUSTCAT fields:

DATESTR:			Date string of earthquake
DATENUM:			Number of days since 1 Jan 0000 (used in matlab)
TYPE:					Event type manually classified by Gary Gibson (e.g. local event, blast, coal blast, teleseismic)
DEPENDENCE:		Describes dependence of event as classified by Gary Gibson (e.g. mainshock, aftershock, foreshock)
LON:					Preferred event lonitude
LAT:					Preferred event latitude
DEP:					Preferred event depth
LOCSRC:				Preferred location source
PREFMW:				Preferred observed MW
PREFMWSRC:		Preferred source of MW
PREFMS:				Preferred observed MS
PREFMSSRC:		Preferred source of MS
PREFmb:				Preferred observed mb
PREFmbSRC:		Preferred source of mb
PREFML:				Preferred observed ML
PREFMLSRC:		Preferred source of ML
REVML:				Revised ML using OBSML corrected using logic described above
OTHERM:				Other magnitude type not captured above
OTHERMTYPE:		Type of other magnitude (e.g., MD, MP)
OTHERMSRC:		Source of other magnitude type
MX_ORIGML:		Takes preferred magnitude type and preserves original ML (PREFML)if ML preferred type.  Field also merges other magnitude types assumed to be equivalent to ML (e.g. MP, M?)
MX_REVML:			Same as MX_ORIGML, but applies ML corrections to original ML (PREFML)
MX_ORIGMLSRC: Source of MX_ORIGML
MS2MW:				MS converted to MW using relations of Scordilis (2006) where applicable
mb2MW:				mb converted to MW using relations of Scordilis (2006) where applicable
ML2MW:				ML converted to MW using relations of Allen (2011) for Australian events (CWA, EA & MLFR)
ML2MWG:				ML converted to MW using relations of Grunthal etal (2009) for all events
PREFMW:				Preferred MW using logic above
PREFMWSRC:		Source of preferred MW
COMM:					Any comments from GG-Cat or ISC associated with event - usually location

******************************************************************************************************************* 
References:

Allen, T. I., T. Dhu, P. R. Cummins, and J. F. Schneider (2006). Empirical attenuation of ground-motion spectral amplitudes in southwestern Western Australia, Bull. Seism. Soc. Am. 96, 572–585.

Allen, T. I., P. R. Cummins, T. Dhu, and J. F. Schneider (2007). Attenuation of ground-motion spectral amplitudes in southeastern Australia, Bull. Seism. Soc. Am. 97, 1279–1292.
