#include "cin.h" 

#pragma warning( push )
#pragma warning( disable : 4042 )		// Disable Warning ('identifier' : has bad storage class)

void LoadBinFile(ErrorCluster *error, const TString &Path);		// Read Binary File in
int32Array ConvertNumber(ErrorCluster *error, uInt8Array Data); // Used to convert Header Info Values
float32Array GetChannelData(ErrorCluster *error, int32 Type, int ChannelIndex, float64 fCutOff);

float32Array LowPass_Filter(ErrorCluster *error, float64 fSample, float64 fCutOff, float32Array InData);
void prewarp(double *a0, double *a1, double *a2, double fc, double fs);
void bilinear(double a0, double a1, double a2, double b0, double b1, double b2, double *k, double fs, float *coef);
void szxform(double *a0, double *a1, double *a2, double *b0, double *b1, double *b2, double fc, double fs, double *k, float *coef);
void bilinear(double a0, double a1, double a2, double b0, double b1, double b2, double *k, double fs, float *coef);

void GetPaceData(ErrorCluster *error, int32 Type, float32 PaceTheshold);
void GetRF_Data(ErrorCluster *error, TStringArray ChannelNames, float64 fCutOff);
float32Array LowPass_Filter(ErrorCluster *error, float64 fSample, float64 fCutOff, float32Array InData);
void GetHeaderInfo(ErrorCluster *error);	// Binary Data Header Information
TString PaceReport(ErrorCluster *error);	// Pace Result out
TString RF_Report(ErrorCluster *error);		// RF_Rectification Result out

	typedef struct {
	int32 dimSizes[2];	// 2-D field sizes
	float32 elt[1];	    // element of the collection
	} TD2;

/* 
 *  INFORMATION STRUCTURES AND VARIABLE DEFINITIONS FOR READ BIN FILE ROUTINES
 *  
 *	Structures used to contain the header info created by: Create Header for Binary Stream to File.vi
 *	For this reason, changes to that VI may impact this type def.
 *
*/

typedef struct {
					LStrHandle channel;			//  Contains the set of analog input channels.
					float32 upperInputLimits;	//  Contains highest expected level of the signals measured.
					float32 lowerInputLimits;	//  Contains lowest expected level of the signals measured.
					float32 range;				//  Onboard input voltage range for a channel.
					uInt16 polarity;			//  Polarity setting for the channel.
					float32 gain;				//  Gain setting for the channel
					uInt16 coupling;			//  Manner in which a signal is connected from one location to another.
					uInt16 inputMode;			//  Input mode for the channel.
					float32 scaleMultiplier;	//  The value used for m in the equation y = mx + b 
												//  where x is the binary reading and y is the resulting scaled value.
					float32 scaleOffset;		//  The value used for b in the equation y = mx + b 
												//  where x is the binary reading and y is the resulting scaled value.
				} GroupCHsettings;
		
typedef struct {
					int32 dimSize;
					GroupCHsettings TypeChannelGroupSetting[50];
				} Channels;
	
Channels DA;
		
int32 HeaderSize;			// Size of the header
int32 FileSize;				// Total File size
float ScanRate;				// scan rate in file
int NumberofChannels = 0;	// Number of channels in file
char *FileBuffer;			// Buffer containing raw (unparsed) lines
char *Rd;					// Used to Parse Filebuffer

// Data for Results
float32Array RV_Pace_Results;	// for XML Results of all RV Pace Data
float32Array A_Pace_Results;	// for XML Results of all Atrial Pace Data
TStringArray RF_CH_Names;		// All rf rectification channel names
float32Array RF_CH_Data;		// All rf rectification Max data

int32Array PaceStart;			// Used to do pace blanking over rf_channels
int32 PaceBlankWidth;			// Data size removed for blanking

TString ReturnString;			// XML Results That will be reurned to user

extern "C" CIN MgErr CINLoad(RsrcFile rf) 
{
    FileBuffer = NULL;
	return noErr;
}

extern "C" CIN MgErr CINUnload(void)
{
	delete [] FileBuffer;
    FileBuffer = NULL;
	return noErr;
}

/// @addtogroup Main
/// @{
/// @addtogroup CINRun Main Function definition
/// @{

/** @brief <B>CINRun:\n</B>
	This is the main entry point into the code.  This code will capture and report RF Rectification\n
	information for the pathways supplied.  This code is intended to be used to read binary file\n
	formats created by 'WriteBinData.cpp'.  The 'Path' input to this function must contain the path\n
	of the binary file being processed.  When 'Type' is set to '0', a Max value for each of the pathways\n
	will be captured.  A lowpass filter may be used on this data.  If a type greater then '0' is selected,\n
	pacing will be looked for and when/if found will be blanked from the RF channels data.  Pacing pulses\n
	width (max/min), amplitude (max/min) and interval (max/min) data will also be returned.  When a pacing\n
	type is selected, the last 2 channels of the binary stream to file data will be considered the pacing\n
	channels.  The second to last channel is the RV and the last is the Atrial.  The following is a more\n
	detailed description of the 'Type' input:\n\n

	Type =\n Return or process type.  Two Types are 'RF' and 'PACE' selected by entering 0-3.\n\n
		   0 = RF \n- This will only return Max waveform data w/o any blanking.  Pace data\n
					returned will contain a result of '-1.0'\n\n
		   1,2,3 = PACE\n - The pace will be detected and processed.  This data must be contained within\n 
						the last two channels of the binary file.  The paces detected will result in\n
						blanking occuring across the rf channels to remove any cross channel issues\n
						do to the pace. A '1' indicates one channel of RV PACE data and a '2' indicates\n
						one channel of Atrial PACE data. '3' indicates both RV and Atrial PACE data.\n\n

	The RF data returned will contain data for each pathway back to case and all other interactions.\n


@param[in,out] error Error Status Information
@param[in] Path	Binary Stream File Path
@param[in] LVChannelNames Array of RF channel names
@param[in] Type Data processing type (Pacing and RF)
@param[in] PaceTheshold Pace detection threshold when type is greater than '0'
@param[in] CutoffFreq digital lowpass filter use on RF channels only. No filter if set less than '0.00'.
@param[out] DataOut				XML formatted data out of all Max RF data and required Pacing information.
*/ 

extern "C" CIN MgErr CINRun(LV_ErrorCluster *error, LStrHandle Path, LV_stringArrayHandle LVChannelNames, int32 *Type, float32 *PaceTheshold, float64 *CutoffFreq, LStrHandle DataOut) 
{

	if (error->status != LVBooleanFalse)	// Check for Error In
		return noErr;	

	ErrorCluster Err(error);					// Declare C++ 'ErrorCluster' Class
	TString File = Path;						// Binary Data File
	TStringArray ChannelNames(LVChannelNames);	// Declare C++ string Array
	ReturnString = TString();				// XML DataOut
	PaceBlankWidth = 0;						// Data size removed for blanking

	//  Load the Binary Data File.  Path is always needed.
	if(File.IsEmpty())
		Err.SetError(-1,"Process RF Data: Need Binary Data Path.");
	else
	{
		LoadBinFile(&Err, File);
		GetHeaderInfo(&Err);
	}

	float32Array TempData = GetChannelData( &Err, 0, 0, -1.00);
	int row = ChannelNames.GetSize();
	int column = TempData.GetSize();
	
	// Check for PACE data because PACE banking will need to occur in RF Data
	RV_Pace_Results.Initialize(-1.0, 6);	// Init RV Data, see PaceReport for format
	A_Pace_Results.Initialize(-1.0, 6);		// Init Atrial Data, see PaceReport for format

	if(!Err.Status() && (*Type < 0 || *Type > 3))
		Err.SetError(-1,"Process RF Data: Invalid Type Range (0 - 3)");
	// There is PACE Data (Note: Pacing does not use a filter)
	else if(!Err.Status() && (*Type > 0))
		GetPaceData(&Err, *Type, *PaceTheshold);
	else
		PaceStart.Initialize(0,0);				// Pace Start of Blanking of RF Data
	
	// Process the RF Data
	if(!Err.Status())
		GetRF_Data(&Err, ChannelNames, *CutoffFreq);

	// Release some memory (so FN LV can consume it). This code doesn't need the file anymore.
	delete [] FileBuffer;
    FileBuffer = NULL;

	// Set/build result out XML. This report is built regardless of an error. 
	ReturnString += PaceReport(&Err);	// Pace Result out
	ReturnString += RF_Report(&Err);	// RF_Rectification Result out
	
	ReturnString.SetLV_String(&Err, &DataOut);	// Update LabVIEW TString
	return Err.LVErr(error);					// Update LabVIEW Error Cluster
}	
/// @}
/// @}

/// @addtogroup File
/// @{
/// @addtogroup LoadBinFile Load Binary data from File.
/// @{

/** @brief <B>LoadBinFile:\n</B>
   This function will load the binary file. This file must be in the format created by
   'WriteBinData.cpp'.  The data in the file will be placed in a buffer to be used by
   other functions.  This buffer will be released when the CIN completes.

@param[in,out] error Error Status Information
@param[in] Path to Binary file to be loaded
*/
void LoadBinFile(ErrorCluster *error, const TString &Path)
{
	if (error->Status())			// if error in, bail out
		return;

	if(FileBuffer != NULL)
	{
		delete [] FileBuffer;
		FileBuffer = NULL;
	}
	
	FILE *in = fopen( Path, "rb" );	// open bin file (read only, binary data)
	if( in == NULL )      
    {
		error->SetError(-1, "LoadBinFile: Unable to open Bin File '%s'", Path.CStr() );
		return;
	}

	if( (FileSize =	_filelength(_fileno(in)) ) <= 0 )	// Get file length
	{
		error->SetError(6, "LoadBinFile: Unable to read file");
		fclose (in);
	}
	else
	{
		if( (FileBuffer = new char[FileSize]) == NULL )	// allocate memory for complete file
			error->SetError( 61, "LoadBinFile: Insufficient memory available");
		else
			fread( FileBuffer, 1, FileSize, in );	// read complete file

		if( ferror(in) != 0 )	// check for file error
		{
			error->SetError(6, "LoadBinFile: Error reading file");
			delete [] FileBuffer;
			FileBuffer = NULL;
		}

		fclose(in);				// Close file
	}
}
/// @}
/// @}

/// @addtogroup File
/// @{
/// @addtogroup GetHeaderInfo Get Binary data header Information.
/// @{

/** @brief <B>GetHeaderInfo:\n</B>
   This function will parse the binary data header information. The information contained
   in this header provides the scan rate used to collect the data, Inter channel delay as well
   as the group channel settings.  The information read by this function was created by the
   'Create Header for Binary Stream to File.vi'.

@param[in,out] error Error Status Information
*/	
void GetHeaderInfo(ErrorCluster *error)
{
	if (error->Status())	// if error in, bail out
		return;

	int32 ChannelListSize;				// Number of Characters in Channel list
	float InterChannelDelay;			// InterChannelDelay
	uInt8Array HInfo;					// Header information
	HInfo.Initialize(0,8);
	Rd = FileBuffer;

	union datum   // Declare a union that can hold the following:
	{
		int32 i;  // int value
		float f;  // float value
	};

	union datum u;

	for(int i = 0; i < 8; i++, Rd++)	// bytes 0-3 are the size of the Header information
		HInfo[i] = *Rd;					// bytes 4-7 are size of channel list
	
	HeaderSize = ConvertNumber(error, HInfo)[0];
	ChannelListSize = ConvertNumber(error, HInfo)[1];
	long ChGroupSettingSize;
	int DataIndex;	

	HInfo.Initialize(0, HeaderSize);				// Re-Init for proper header size
	
	for(i = 0; i < ChannelListSize; i++, Rd++);		// jump over channel list

	for(i = 0; i < 4; i++, Rd++)	// Channel group setting size is the next 4 bytes
		HInfo[i] = *Rd;

	ChGroupSettingSize = ConvertNumber(error, HInfo)[0];

	for(i = 0; i < 4; i++, Rd++)	// Get number of channels
		HInfo[i] = *Rd;

	NumberofChannels = ConvertNumber(error, HInfo.SubArray(0, 4))[0];
			
	int32Array GroupSettings;
	
	for( int z = 0; z < NumberofChannels; z++, Rd++)
	{		
		for(i = 0; i < 4; i++, Rd++)	// Get new Dataindex
			HInfo[i] = *Rd;

		DataIndex = ConvertNumber(error, HInfo.SubArray(0, 4))[0];
		
		for(i = 0; i < DataIndex; i++, Rd++);	// Jump Text String

		for(i = 0; i < 40; i++, Rd++)			// Get channel group setting 
			HInfo[i] = *Rd;

		GroupSettings = ConvertNumber(error, HInfo);
		
		u.i = GroupSettings[0];
		DA.TypeChannelGroupSetting[z].upperInputLimits = u.f;
		u.i = GroupSettings[1];
		DA.TypeChannelGroupSetting[z].lowerInputLimits = u.f;
		u.i = GroupSettings[2];
		DA.TypeChannelGroupSetting[z].range = u.f;
		DA.TypeChannelGroupSetting[z].polarity = (uInt16)(GroupSettings[3] >> 16 );
		u.i = (GroupSettings[3] & 0xFFFF) << 16 | ((GroupSettings[4]) >> 16 & 0xFFFF);
		DA.TypeChannelGroupSetting[z].gain = u.f;
		DA.TypeChannelGroupSetting[z].coupling = (uInt16)(GroupSettings[4] & 0xFFFF);
		DA.TypeChannelGroupSetting[z].inputMode = (uInt16)(GroupSettings[5] >> 16 );
		u.i = (GroupSettings[5] & 0xFFFF) << 16 | ((GroupSettings[6]) >> 16 & 0xFFFF);
		DA.TypeChannelGroupSetting[z].scaleMultiplier = u.f;	
		u.i = (GroupSettings[6] & 0xFFFF) << 16 | ((GroupSettings[7]) >> 16 & 0xFFFF);	
		DA.TypeChannelGroupSetting[z].scaleOffset = u.f;
		/*  scan rate in file is not always valid */		
		u.i = (GroupSettings[7] & 0xFFFF) << 16 | ((GroupSettings[8]) >> 16 & 0xFFFF);
		ScanRate = u.f;
		/*  InterChannelDelay is not always valid */
		u.i = (GroupSettings[8] & 0xFFFF) << 16 | ((GroupSettings[9]) >> 16 & 0xFFFF);
		InterChannelDelay = u.f;													
		
		Rd = Rd - 11;	//  ScanRate & InterChannelDelay are only valid for the last Group Setting
	}
}
/// @}
/// @}

/// @addtogroup ProcessChannelData
/// @{
/// @addtogroup GetChannelData Get channel data.
/// @{

/** @brief <B>GetChannelData:\n</B>
  This is the DAQ channel data.  This function converts the Binary data (format of LITTLE ENDIAN)\n
  to a float32 Array.  This function will also filter the data if requested.  The returned data\n
  will remove the first 0.01 seconds of the array.  This 'time' that is removed is data required\n
  to create/start the history for the filter.  Since a Lowpass digital filter is expected to be\n
  used to process rectification, all data must include this 10 msec addition (2 sec wanted = 2.01).\n
  The pace channels must do the same thing even when filtering isn't done.  This is because blanking\n
  across the rf channels will occur and aligning the index for blanking is easier if you don't\n
  have to account for variations in the start and stop.  The returned data will also have the\n
  pace blanking data removed from the rf data.  This blanking occurs after the entire DAQ capture\n
  has been run through the filter.  The blanking is deleted from the array.  This is to speed up\n
  the cross channel min and max calcs (smaller arrays = fewer subtractions to make = faster execution).\n\n   

@param[in,out] error Error Status Information
@param[in] Type Data processing type (1-3 Pacing and 0 RF)
@param[in] ChannelIndex DAQ channel index
@param[in] fCutOff digital lowpass filter use on RF channels only. No filter if set less than '0.00'.
*/	
float32Array GetChannelData( ErrorCluster *error, int32 Type, int ChannelIndex, float64 fCutOff)
{		
		if (error->Status())	// if error in, bail out
			return float32Array();

		float32Array TempData;	// Temp Data for filter and pre-trigger removal

		//  FYI:  This function is called and controlled internally so some of the error
		//  handling that was in the daq code is not needed and would simply add to the
		//  execution time.

		uInt8Array charData;
		charData.Initialize(0, (FileSize - HeaderSize - 4) / NumberofChannels);

		Rd = FileBuffer;
		for(int p =0; p < HeaderSize + 4 ; p++, Rd++);	//  find data remember first 4 bytes are the header size
		
		int a = 0;
		if(ChannelIndex != 0)
			for(a = 0; a < ChannelIndex * 2; a++,Rd++);	//  Find channel index

		for(int i = 0; i < charData.GetSize(); i++)
		{
			charData[i] = *Rd;
			i++;
			Rd++;
			charData[i] = *Rd;
			for(a = 0; a < NumberofChannels * 2 - 1; a++,Rd++);
		}

		int16Array values;
		values.Initialize(0,charData.GetSize()/2);
		
		// Get binary data (LITTLE ENDIAN format)
		for(i = 0, a = 0; i < values.GetSize(); i++, a +=2)
			values[i] = (int16)((uInt16)charData[a+1] << 8 | charData[a]);
			
		float32Array VData;
		VData.Initialize(0.0, values.GetSize());

		// Convert the binary data to volts (binary file value to float32 array volts).  
		// All of this data should be bi-polar for this testing.
		if( DA.TypeChannelGroupSetting[0].polarity == 1 /*&& Type == 0*/) // RF Channel or Pace done the same 
		{
			for( i = 0; i < VData.GetSize(); i++)
			{
				// bi-polar
				VData[i] = (float32)(values[i] * DA.TypeChannelGroupSetting[ChannelIndex].scaleMultiplier + 
				DA.TypeChannelGroupSetting[ChannelIndex].scaleOffset);
			}
		}
		else
		{
			error->SetError(-1, "Process RF Data: uni-polar data not accepted. Check DAQ config.");
			return float32Array();
		}
			
		// Make sure there is a reasonable amount of data
		int StartIndex = int(ScanRate * 0.01);  // 10 msec pre-trigger

		if( VData.GetSize() <= StartIndex)
		{
			error->SetError(-1, "Process RF Data: Not enough data. Need > 10 msec.");
			return float32Array();
		}

		// Now lets filter the data if requested.
		if( fCutOff > 0.00)
		{
			float32Array TData = LowPass_Filter(error, ScanRate, fCutOff, VData);
			
			if (error->Status())	// if error, bail out
				return float32Array();
			else
				TempData = TData.SubArray(StartIndex, -1);	// filtered data - pre-trigger
		}
		else
			TempData = VData.SubArray(StartIndex, -1);		// no filter data - pre-trigger

		// Lets do the blanking
		if((Type == 0) && (PaceBlankWidth > 0))  // don't need to do this if it is a pace or no blanking
		{
			// There is only one Pace Blanking Width.  It is the largest value BlankingWidth captured.
			// In order to capture the cross channel min and max values, the two channels (arrays) being 
			// compared will be subtracted from one another.  To reduce the amount of time this takes,
			// the banking will be removed from the arrays.
						
			// There may be a partial pace at the beginning or end that need to be blanked.
			int ArrayStart = PaceStart[0]/2;
			int ArrayEnd = (TempData.GetSize() - PaceStart[PaceStart.GetSize()])/2;

			// First intialize the new array
			float32Array BlankedData;
			// The size is the current array - a pace blanking width for each pace detected
			int BlankedDataSize = TempData.GetSize() - ArrayStart - ArrayEnd - (PaceStart.GetSize() * PaceBlankWidth);
			BlankedData.Initialize(0.0, BlankedDataSize);

			// TempData will control the for loop but we will also need to control the PaceStart
			// index and the index we are placing data into the new BlankData Array
			int PaceStartIndex = 0;
			int BlankDataIndex = 0;
			bool getenddata = false; // Protect PaceStartIndex for over indexing array
			
			// Here is the build of the blank array.  Both cases (if/else) are controlled
			// to prevent an invalid array index issue from occuring
			for(i = ArrayStart; (i < TempData.GetSize()) && (BlankDataIndex < BlankedDataSize); )
			{
				if(!getenddata)
				{
					if( i == PaceStart[PaceStartIndex]) // This is blanked
					{	
						i += PaceBlankWidth;	// Index TempData
						PaceStartIndex++;		// Get Next Blank location
						if(PaceStartIndex >= PaceStart.GetSize()) 
							getenddata = true;	// Ran out of Paces
					}
					else
					{
						BlankedData[BlankDataIndex] = TempData[i];	// Keep Data
						BlankDataIndex++;	// Index to next storage location
						i++;				// Index next data
					}
				}
				else
				{
					BlankedData[BlankDataIndex] = TempData[i];	// Keep Data
					BlankDataIndex++;	// Index to next storage location
					i++;				// Index next data
				}

			}
			return BlankedData; // return the blanked data
		}
		else
			return TempData;	// no blanking, so return
}
/// @}
/// @}

/// @addtogroup Utility
/// @{
/// @addtogroup ConvertNumber Convert 8 bit to 32.
/// @{

/** @brief <B>ConvertNumber:\n</B>
   This function converts the LV header numerical values to 
   an Int32.  The format being converted from that file type
   is BIG ENDIAN.
   
@param[in,out] error Error Status Information
@param[in] Data 8-bit array to be converted
*/
int32Array ConvertNumber(ErrorCluster *error, uInt8Array Data)
{
	if (error->Status())  // if error in, bail out
		return int32Array();

	if(Data.GetSize() < 4)
	{
		error->SetError(-1, "ConvertNumber: No data to convert.");
		return int32Array();
	}
	
	int32Array ConvData;
	ConvData.Initialize(0, Data.GetSize() / 4);

	// BIG ENDIAN
	for(int i = 0, a = 0; i < Data.GetSize() / 4; i++, a = a + 4) 
		ConvData[i] = (int32)((int32)(Data[a] << 24) | (int32)(Data[a+1] << 16) | (int32)(Data[a+2] << 8) | (int32)(Data[a+3]));
	
	return ConvData;
}
/// @}
/// @}

/// @addtogroup ProcessChannelData
/// @{
/// @addtogroup GetPaceData Get Pace data.
/// @{

/** @brief <B>GetPaceData:\n</B>
 This gets all the pacing data.  It needs to be called prior to calling GetRF_Data() because paces\n 
 will be blanked in the RF data.  Two Pacing result arrays can be populated by this function.\n 
 RV_Pace_Results will be created if Type is a 1 or 3.	A_Pace_Results will be created if Type is\n
 a 2 or 3.  The pace results for each array are defined as follows:\n\n
  Pace_Results[0] = PW_Max\n
  Pace_Results[1] = PW_Min\n   
  Pace_Results[2] = Ampl_Max\n  
  Pace_Results[3] = Ampl_Min\n  
  Pace_Results[4] = Interval_Max\n  
  Pace_Results[5] = Interval_Min\n\n 
 
 The other thing this function does is to capture an array of blanking with respect to the pace.\n
 'PaceStart' start array contains the locations (by index) of where the pace banking should occur.\n
 The start of blanking is 1 mSec prior to the pace edge being detected.  The blanking is calculated\n
 as follows:\n\n 
 PaceBlankWidth = 1 mSec + Pace PW + 6.8 mSec + 1 mSec \n
 Where:\n
 1 mSec = 1 mSec prior to the pace edge and 1 mSec after the end of recharge.\n
 6.8 mSec = nominal recharge + interchannel delay.\n\n
 
 The blanking used will be with respect to the largest PW measured. It probably should have\n
 been an array for each channel but I didn't think about it until later.  It is hard to say how much\n
 it really matters.\n\n 

 Note: Code Speed is a product of instruction cycles executed to perform the function not code size.\n\n 

@param[in,out] error Error Status Information
@param[in] Type Data processing type (Pacing and RF)
@param[in] PaceTheshold Pace detection threshold when type is greater than '0'
*/	
void GetPaceData(ErrorCluster *error, int32 Type, float32 PaceTheshold)
{
	if (error->Status())  // if error in, bail out
		return;

	// Going to need a couple of temp arrays to get some of the data figured out
	int32Array Temp_VP_Result;		// RV
	Temp_VP_Result.Initialize(0,0);
	int32Array Temp_AP_Result;		// Atrial
	Temp_AP_Result.Initialize(0,0);

	// Initialize items for pace calcs
	float32 PW_max = -1.0;
	float32 PW_min = -1.0;
	float32 Ampl_max = -1.0;
	float32 Ampl_min = -1.0;
	float32 Int_max = -1.0;
	float32 Int_min = -1.0;

	int32Array Start_Array_RV;	// PaceStart for RV		
	int32Array Start_Array_A;	// PaceStart for the Atrial

	int32Array TempPaceData;

	int Index = 0; // Channel Index in Binary File

	// Do a width/threshold search for pace edges
	int ThresholdWidth = int(ScanRate * 0.00005);	// 50 samples w/1M sample rate
	int WidthCount = 0;								// Used to track samples

	// Values needed for PaceBlankWidth = 1 mSec + Pace PW + 6.8 mSec + 1 mSec
	int PaceBlankExtra = int(ScanRate * 0.001);		// 1 mSec per BlankWidth Calc
	int RechargeWidth = int(ScanRate * 0.0068);		// 6.8 mSec per BlankWidth Calc

	// Value used to calculate peak pace amplitude
	int ten_us = int(ScanRate * 0.00001);			// 10 uSec	

	// Type = 1 is Vent (RV) Data only.  Type = 3 is RV and Atrial Both
	if(Type == 1 || Type == 3)
	{
		Index = NumberofChannels - 2;	// RV (2nd to last) Channel in Binary File

		// Get RV Data, Any PACE Data DAQ should capture 2 channels. Type
		// is used to determine which of the 2 are processed.

		// Get the data, No filter is used on Pace Data
		float32Array Temp = GetChannelData( error, Type, Index, -1.0);

		if (error->Status())  // if error, bail out
			return;

		int TempSize = Temp.GetSize();

		// Need to capture pace detected start and stop edges.  These
		// arrays will have the data Appended which is slow but most
		// testing will not have more than 10 paces captured.
		int32Array Temp_Start;
		int32Array Temp_End;

		// Find the Paces.  This could be a function call but saves instruction
		// cycles leaving it here.  It toggles between looking for the falling
		// edge and the rising edge of the pace.  The pace should be negative.
		for( int i = 0; i < TempSize; i++)
		{
			// Do a threshold search. We are looking for the negative leading
			// pace edge.  Because of noise concerns, we want a number of
			// points below this threshold.
			if( PaceTheshold >= Temp[i])
				WidthCount++;
			else
				WidthCount = 0;

			if( WidthCount >= ThresholdWidth)
			{
				Temp_Start.Insert((i - ThresholdWidth), -1); // Get Index
				WidthCount = 0;	// Reset Count
				
				// Now look for rising edge of Pace
				for( ;(WidthCount < ThresholdWidth) && (i < TempSize); i++)
				{
					if( PaceTheshold <= Temp[i])
						WidthCount++;
					else
						WidthCount = 0;

					if( WidthCount == ThresholdWidth)
						Temp_End.Insert( i - ThresholdWidth, -1); // Get Index
				
				}
				
				// FYI: The first if/else should reset WidthCount
			}
			
		}

		WidthCount = 0;	// Reset Count 

		//  If there is data then calcs may be done
		if(Temp_Start.GetSize() > 0 && Temp_End.GetSize() > 0)
		{
			// Initialize items for pace calcs
			PW_max = -1.0;
			PW_min = -1.0;
			Ampl_max = -1.0;
			Ampl_min = -1.0;
			Int_max = -1.0;
			Int_min = -1.0;

			// Initialize storage arrays
			float32Array PW_Array;
			PW_Array.Initialize( -1.0, 0);		// PW
			float32Array Ampl_Array;
			Ampl_Array.Initialize( -1.0, 0);	// Amplitude
			float32Array Int_Array;
			Int_Array.Initialize( -1.0, 0);		// Interval
			
			int End = 0;	// Used to set pace array size End
			int Start = 0;	// Index first used to set Pace Array Start

			// Initialize End to the shorter of two arrays
			if( Temp_Start.GetSize() <= Temp_End.GetSize())
				End = Temp_Start.GetSize();
			else
				End = Temp_End.GetSize();

			// Make sure there is enough pre-pace data
			if( (Temp_Start[0] - PaceBlankExtra) < 0)
				Start = 1;	// Skip leading and trailing pace
			// Make sure there is enough post-pace data 
			if( (Temp_End[TempPaceData.GetSize()] + PaceBlankExtra) > Temp.GetSize())
				End--;
			
			int sizes = End - Start; // Size of PW, Ampl, Interval arrays

			if(sizes > 0) // There is at least one pace
			{
				// We can do at least some pacing measurements
				PW_Array.Initialize( -1.0, sizes);		// PW
				Ampl_Array.Initialize( -1.0, sizes);	// Amplitude
				Start_Array_RV.Initialize( -1.0, Temp_Start.GetSize());// Will have -1 mSec

				// Need more then one for an Interval
				if(sizes > 1)
				{
					// For 6 paces there are only 5 intervals
					Int_Array.Initialize( -1, sizes);	// Interval
					int b = 0;

					for( i = Start; i < End; i+=2, b++)
						Int_Array[b] = (float32)(Temp_Start[i+1] - Temp_Start[i]);
				}
				// Get the PW and Pace blanking start locations in RV Pace Data
				int b = 0;
				for( i = Start; i < End; i++, b++)
					PW_Array[b] = (float32)(Temp_End[i] - Temp_Start[i]);	// PW
				
				// used for blanking
				int Index = 0;
				for(i = 0; i < Temp_Start.GetSize(); i++)
				{
					Index = Temp_Start[i] - PaceBlankExtra;
					if(Index < 0)
						Start_Array_RV[i] = 0;
					else
						Start_Array_RV[i] = Index;		
				}
				
				// Get Peak Pace Amplitudes.  This is captured by taking a 10 uSec mean
				// of the data located 10 uSec in from the pace edge
				b = 0;
				for( i = Start; i < End; i++, b++)
				{
					float32 temp_pk = 0.0;					// Used to Calculate Peak
					int pk_s = Temp_Start[i] + ten_us;		// Peak pace start (10 us from edge)
					int pk_e = pk_s + ten_us;				// Peak pace End ( 10 us average)

					// Find peak average
					for( int a = pk_s; a < pk_e; a++)
						temp_pk += Temp[a];
					
					Ampl_Array[b] = (float32)(temp_pk/(float32)ten_us);// Capture average
				}

				// Get Max, Min, Mean.  These arrays should be typically less then 10
				// values so I'm going to write this in an easy way to follow rather
				// then 'fast' execution.  Program will be large but it isn't a concern.
				for(i = 0; i < sizes; i++)		// all Arrays are = size
				{

					// Look for max and min value
					if( i == 0)
					{
						PW_max = PW_Array[i];		// Reset PW_max 
						PW_min = PW_Array[i];		// Reset PW_min
						Ampl_max = Ampl_Array[i];	// Reset Ampl_max 
						Ampl_min = Ampl_Array[i];	// Reset Ampl_min
						Int_max = Int_Array[i];		// Reset Int_max 
						Int_min = Int_Array[i];		// Reset Int_min
					}
					else
					{
						if(PW_Array[i] > PW_max)	// Look for PW_Max
							PW_max = PW_Array[i];
						if(PW_Array[i] < PW_min)	// Look for PW_Min
							PW_min = PW_Array[i];

						if(Ampl_Array[i] < Ampl_max)// Look for Ampl_Max (more neg)
							Ampl_max = Ampl_Array[i];
						if(Ampl_Array > Ampl_min)	// Look for Ampl_Min (more pos)
							Ampl_min = Ampl_Array[i];

						if(Int_Array[i] > Int_max)	// Look for Int_Max
							Int_max = Int_Array[i];
						if(Int_Array < Int_min)		// Look for Int_Min
							Int_min = Int_Array[i];
					}
				}
				// Store the data as follows:
				//  Pace_Results[0] = PW_Max
				//  Pace_Results[1] = PW_Min   
				//  Pace_Results[2] = Ampl_Max  
				//  Pace_Results[3] = Ampl_Min   
				//  Pace_Results[4] = Interval_Max  
				//  Pace_Results[5] = Interval_Min  
				RV_Pace_Results[0] = (PW_max * 1000.0)/ScanRate;
				RV_Pace_Results[1] = (PW_min * 1000.0)/ScanRate;
				RV_Pace_Results[2] = Ampl_max;
				RV_Pace_Results[3] = Ampl_min;
				RV_Pace_Results[4] = (Int_max/ScanRate) * 1000.0;
				RV_Pace_Results[5] = (Int_min/ScanRate) * 1000.0;

				// Set the Pace Blank width.  Let the Atrial Pacing section
				// figure out if a larger value should be used.  This code
				// should be changed to better handle different PWs but the
				// testing will be done with 1 and 1.5 mSec timing.  If one
				// Pace channel has an additional 500 us of data blanked it
				// isn't going to matter much.
				// PaceBlankWidth = 1 mSec + Pace PW + 6.8 mSec + 1 mSec
				PaceBlankWidth = PaceBlankExtra+PW_max+RechargeWidth+PaceBlankExtra;
			}
		}		
	}
			
	// Type = 2 is Atrial Data only.  Type = 3 is RV and Atrial Both
	if(Type == 2 || Type == 3)
	{
		Index = NumberofChannels - 1;	// Atrial (last) Channel in Binary File
		
		// Get Atrial Data, Any PACE Data DAQ should capture 2 channels. Type
		// is used to determine which of the 2 are processed.

		// Get the data, No filter is used on Pace Data
		float32Array Temp = GetChannelData( error, Type, Index, -1.0);

		if (error->Status())  // if error, bail out
			return;

		int TempSize = Temp.GetSize();

		// Need to capture pace detected start and stop edges.  These
		// arrays will have the data Appended which is slow but most
		// testing will not have more than 10 paces captured.
		int32Array Temp_Start;
		int32Array Temp_End;

		// Find the Paces.  This could be a function call but saves instruction
		// cycles leaving it here.  It toggles between looking for the falling
		// edge and the rising edge of the pace.  The pace should be negative.
		for( int i = 0; i < TempSize; i++)
		{
			// Do a threshold search. We are looking for the negative leading
			// pace edge.  Because of noise concerns, we want a number of
			// points below this threshold.
			if( PaceTheshold >= Temp[i])
				WidthCount++;
			else
				WidthCount = 0;

			if( WidthCount >= ThresholdWidth)
			{
				Temp_Start.Insert((i - ThresholdWidth), -1); // Get Index
				WidthCount = 0;	// Reset Count
				
				// Now look for rising edge of Pace
				for( ;(WidthCount < ThresholdWidth) && (i < Temp.GetSize()); i++)
				{
					if( PaceTheshold <= Temp[i])
						WidthCount++;
					else
						WidthCount = 0;

					if( WidthCount == ThresholdWidth)
						Temp_End.Insert( i - ThresholdWidth, -1); // Get Index
				
				}
				
				// FYI: The first if/else should reset WidthCount
			}
			
		}

		WidthCount = 0;	// Reset Count

		//  If there is data then calcs may be done
		if(Temp_Start.GetSize() > 0 && Temp_End.GetSize() > 0)
		{
			// Initialize items for pace calcs
			PW_max = -1.0;
			PW_min = -1.0;
			Ampl_max = -1.0;
			Ampl_min = -1.0;
			Int_max = -1.0;
			Int_min = -1.0;

			// Initialize storage arrays
			float32Array PW_Array;
			PW_Array.Initialize( -1.0, 0);		// PW
			float32Array Ampl_Array;
			Ampl_Array.Initialize( -1.0, 0);	// Amplitude
			float32Array Int_Array;
			Int_Array.Initialize( -1.0, 0);		// Interval
					
			int End = 0;	// Used to set pace array size End
			int Start = 0;	// Index first used to set Pace Array Start

			// Initialize End to the shorter of two arrays
			if( Temp_Start.GetSize() <= Temp_End.GetSize())
				End = Temp_Start.GetSize();
			else
				End = Temp_End.GetSize();

			// Make sure there is enough pre-pace data
			if( (Temp_Start[0] - PaceBlankExtra) < 0)
				Start = 1;	// Skip leading and trailing pace
			// Make sure there is enough post-pace data 
			if( (Temp_End[TempPaceData.GetSize()] + PaceBlankExtra) > Temp.GetSize())
				End--;
			
			int sizes = End - Start; // Size of PW, Ampl, Interval arrays

			if(sizes > 0) // There is at least one pace
			{
				// We can do at least some pacing measurements
				PW_Array.Initialize( -1.0, sizes);		// PW
				Ampl_Array.Initialize( -1.0, sizes);	// Amplitude
				Start_Array_A.Initialize( -1.0, Temp_Start.GetSize());	// Will have -1 mSec

				// Need more then one for an Interval
				if(sizes > 1)
				{
					Int_Array.Initialize( -1, sizes);	// Interval
					int b = 0;

					for( i = Start; i < End; i+=2, b++)
						Int_Array[b] = (float32)(Temp_Start[i+1] - Temp_Start[i]);
				}
				// Get the PW and Pace blanking start locations in RV Pace Data
				int b = 0;
				for( i = Start; i < End; i++, b++)
					PW_Array[b] = (float32)(Temp_End[i] - Temp_Start[i]);	// PW

				// used for blanking
				int Index = 0;
				for(i = 0; i < Temp_Start.GetSize(); i++)
				{
					Index = Temp_Start[i] - PaceBlankExtra;
					if(Index < 0)
						Start_Array_A[i] = 0;
					else
						Start_Array_A[i] = Index;		
				}
				
				// Get Peak Pace Amplitudes.  This is captured by taking a 10 uSec mean
				// of the data located 10 uSec in from the pace edge
				b = 0;
				for( i = Start; i < End; i++, b++)
				{
					float32 temp_pk = 0.0;					// Used to Calculate Peak
					int pk_s = Temp_Start[i] + ten_us;		// Peak pace start (10 us from edge)
					int pk_e = pk_s + ten_us;				// Peak pace End ( 10 us average)

					// Find peak average
					for( int a = pk_s; a < pk_e; a++)
						temp_pk += Temp[a];
					
					Ampl_Array[b] = (float32)(temp_pk/(float32)ten_us);// Capture average
				}

				// Get Max and Min.  These arrays should be typically less then 10
				// values so I'm going to write this in an easy way to follow rather
				// then 'fast' execution.  Program will be large but it isn't a concern.
				for(i = 0; i < sizes; i++)		// all Arrays are = size
				{
					// Look for max and min value
					if( i == 0)
					{
						PW_max = PW_Array[i];		// Reset PW_max 
						PW_min = PW_Array[i];		// Reset PW_min
						Ampl_max = Ampl_Array[i];	// Reset Ampl_max 
						Ampl_min = Ampl_Array[i];	// Reset Ampl_min
						Int_max = Int_Array[i];		// Reset Int_max 
						Int_min = Int_Array[i];		// Reset Int_min
					}
					else
					{
						if(PW_Array[i] > PW_max)	// Look for PW_Max
							PW_max = PW_Array[i];
						if(PW_Array[i] < PW_min)	// Look for PW_Min
							PW_min = PW_Array[i];

						if(Ampl_Array[i] < Ampl_max)// Look for Ampl_Max (more neg)
							Ampl_max = Ampl_Array[i];
						if(Ampl_Array > Ampl_min)	// Look for Ampl_Min (more pos)
							Ampl_min = Ampl_Array[i];

						if(Int_Array[i] > Int_max)	// Look for Int_Max
							Int_max = Int_Array[i];
						if(Int_Array < Int_min)		// Look for Int_Min
							Int_min = Int_Array[i];
					}
				}
				// Store the data as follows:
				//  Pace_Results[0] = PW_Max
				//  Pace_Results[1] = PW_Min   
				//  Pace_Results[2] = Ampl_Max  
				//  Pace_Results[3] = Ampl_Min  
				//  Pace_Results[4] = Interval_Max  
				//  Pace_Results[5] = Interval_Min  
				A_Pace_Results[0] = (PW_max * 1000.0)/ScanRate;
				A_Pace_Results[1] = (PW_min * 1000.0)/ScanRate;
				A_Pace_Results[2] = Ampl_max;
				A_Pace_Results[3] = Ampl_min;
				A_Pace_Results[4] = (Int_max/ScanRate) * 1000.0;
				A_Pace_Results[5] = (Int_min/ScanRate) * 1000.0;

				// Set the Pace Blank width.  Check to see if the Atrial has
				// a larger width then the Vent.  This value was Intialized to
				// 0 so if the Vent wasn't captured it wont matter. This code
				// should be changed to better handle different PWs but the
				// testing will be done with 1 and 1.5 mSec timing.  If one
				// Pace channel has an additional 500 us of data blanked it
				// isn't going to matter much.
				// PaceBlankWidth = 1 mSec + Pace PW + 6.8 mSec + 1 mSec
				int32 A_PaceBlankWidth = PaceBlankExtra+PW_max+RechargeWidth+PaceBlankExtra;

				if(A_PaceBlankWidth > PaceBlankWidth)
					PaceBlankWidth = A_PaceBlankWidth;
			}
		}
	}

	// Take care of PaceStart (blanking array)
	if( Start_Array_RV.GetSize() == 0 && Start_Array_A.GetSize() == 0)		// No Pace Data
		PaceStart.Initialize(0,0);	// Pace Start of Blanking of RF Data
	else if(Start_Array_RV.GetSize() > 0 && Start_Array_A.GetSize() == 0)	// RV only Data
		PaceStart = Start_Array_RV;
	else if(Start_Array_RV.GetSize() == 0 && Start_Array_A.GetSize() > 0)	// Atrial only Data
		PaceStart = Start_Array_A;
	else	// Have both
	{
		int Size = Start_Array_A.GetSize() + Start_Array_RV.GetSize();
		int b = 0;

		PaceStart.Initialize(0,Size);

		// This is all captured by a software start so the first pace detected will
		// state the order of the array
		if(Start_Array_A[0] > Start_Array_RV[0])	// Then RV is First
		{
			for(int i = 0; i < Size; i+=2, b++)
			{
				if(b < Start_Array_RV.GetSize())
					PaceStart[i] = Start_Array_RV[b];
				if(b < Start_Array_A.GetSize())
					PaceStart[i+1] = Start_Array_A[b];
			}
		}
		else	// Atrial is first
		{
			for(int i = 0; i < Size; i+=2, b++)
			{
				if(b < Start_Array_A.GetSize())
					PaceStart[i] = Start_Array_A[b];
				if(b < Start_Array_RV.GetSize())
					PaceStart[i+1] = Start_Array_RV[b];
			}
		}
	
	}
	
}
/// @}
/// @}

/// @addtogroup ProcessChannelData
/// @{
/// @addtogroup GetRF_Data Get RF data.
/// @{

/** @brief <B>GetRF_Data:\n</B>
 This gets the rf_rectification data.  The data it gets is entirely based upon the\n
 size of the ChannelNames array.  Obviously, the single channel data is a portion\n
 of the data captured.  This is the channel name with respect to case or the reference\n
 The 'respect to case' data is then compared to each of the other paths.  This results\n
 in the following for 2 paths of rvtip and rvring:\n\n
	1. rvtip_to_case\n
	2. rvring_to_case\n
	3. rvtip_to_rvring\n\n
 This function will fill two arrays.  One will contain all the path names in the order\n
 they are processed.  The other will contain an array of the Max followed by Min values.\n
 The Max and Min values will be with respect to zero.  In other words, if a negative\n
 value is further from zero then a positive value it will be reported as a Max.  No\n
 negative values will be returned.  Rather the absolute value.\n\n 

@param[in,out] error Error Status Information
@param[in] ChannelNames Array of RF channel names
@param[in] fCutOff digital lowpass filter use on RF channels only. No filter if set less than '0.00'.
*/	
void GetRF_Data(ErrorCluster *error, TStringArray ChannelNames, float64 fCutOff)
{
	if (error->Status())  // if error in, bail out
		return;

	int RF_CH_Names_Size = 0;	// Get Size of all ChannelName combinations

	// The name size is a summation of all element units. Example:
	//	ChannelNames Array size = 6
	//  RF_CH_Names_Size = 6 + 5 + 4 + 3 + 2 + 1 =  21
	for(int i = ChannelNames.GetSize(); i > 0; i--)
		RF_CH_Names_Size += i;

	// Initialize the RF_CH_Names array (used to build report XML)
	RF_CH_Names.Initialize(TString(), RF_CH_Names_Size);

	// Initialize the RF_CH_Data array (used to build report XML)
	// This data array contains a Max value for each path.
	RF_CH_Data.Initialize(0, RF_CH_Names_Size);

	// The data will be processed using a 2D Array.  Coulmns of the array will
	// be subtracted from each other. All coulmns need to be compared. Here is
	// what it should look like for ChannelNames size of 3.  There will be
	// 6 RF_CH_Data elements (3+2+1 = 6) so the column arrays should index as
	// follows:
	//   A   B
	//  [0] [0] - path to case (don't subtract)
	//  [0] [1]
	//  [0] [2]
	//  [1] [1] - path to case (don't subtract)
	//  [1] [2]
	//  [2] [2] - path to case (don't subtract)
	// When the index value of A and B are equal the Max and Min of the data
	// array should be taken without doing any subtraction.  This data indicates
	// a path to case.

	int32Array Index_A;
	Index_A.Initialize(0, RF_CH_Names_Size); // First Column Index
	int32Array Index_B;
	Index_B.Initialize(0, RF_CH_Names_Size); // Second Column Index

	for( i = 0; i < RF_CH_Data.GetSize(); i++)	// control index storing data
	{
		for( int a = 0; a < ChannelNames.GetSize(); a++)	// control 'A'
		{
			int b = a;										// walk 'B' w/'A'
			for(; b < ChannelNames.GetSize(); b++, i++)		// control 'B'
			{
				Index_A[i] = a;	// Store data
				Index_B[i] = b; // Store data
			}
		}

	}
	
	// Now get the data.  Channel '0' is always used. This can only be
	// RF data.  A filter may be used.  This first call is used to 
	// define a 2D array size.  The format of the array will be each
	// column will contain a different channels data.  The row will
	// be indexed to process the channels array: TWO_D_ARRAY[ROW][COLUMN]
	float32Array TempData = GetChannelData( error, 0, 0, fCutOff);
	int row = TempData.GetSize();
	int column = ChannelNames.GetSize();

	// This uses a sloppy loop method to initialize a 2D array
	float32 **TWO_D_ARRAY = new float32*[column];
	for(i = 0; i < column; ++i) 
		TWO_D_ARRAY[i] = new float32[row];

		// Capture Channel '0' (so you don't have to get the data again)
	for( i = 0; i < TempData.GetSize(); i++)
		TWO_D_ARRAY[0][i] = TempData[i];

	// Now get the rest of the data
	for( int a = 1; a < ChannelNames.GetSize(); a++)	// Increment channel
	{
		TempData = GetChannelData( error, 0, a, fCutOff);
		for( i = 0; i < TempData.GetSize(); i++)
			TWO_D_ARRAY[a][i] = TempData[i];
	}

	// Now get all the Max and Min Values
	int r = 0;	// Row Index
	int c = 0;	// Index for RF_CH_Data (Max data)
	float32 max = 0.0;	// use to find max, 0 won't be used in case of channel offset
	float32 min = 0.0;	// use to find min, 0 won't be used in case of channel offset
	
	// We can also build the RF_CH_Names Array.  The format is ['Path A'_to_'Path B'].
	// When Index 'A' = 'B', the format is ['Path A'_to_Case].
	TString FormatResult;
	TString TO = TString("_to_");
	TString CASE = TString("Case");

	for( i = 0; i < RF_CH_Names_Size; i++, c++)
	{
		a = Index_A[i];
		int b = Index_B[i];

		FormatResult = ChannelNames[a]; // Always use Item in first Array
		FormatResult += TO;				// Add '_to_'

		if( a!=b)	// Not a Path to Case 
		{
			FormatResult += ChannelNames[b];// Add second path
			r = 0;							// reset row index

			TempData[r] = TWO_D_ARRAY[a][r] - TWO_D_ARRAY[b][r];
			max = TempData[r];	// Reset max 
			min = TempData[r];	// Reset min
			
			for(r = 1; r < row; r++)
			{
				TempData[r] = TWO_D_ARRAY[a][r] - TWO_D_ARRAY[b][r];

				if(TempData[r] > max)	// Look for Max
					max = TempData[r];
				if(TempData[r] < min)	// Look for Min
					min = TempData[r];
			}

		}
		else		// This is Path to Case
		{
			FormatResult += CASE;		// Add 'Case'
			r = 0;						// reset row index
			max = TWO_D_ARRAY[b][r];	// Reset max
			min = TWO_D_ARRAY[b][r];	// Reset min

			for(r = 1; r < row; r++)
			{
				if(max < TWO_D_ARRAY[b][r])	// Look for Max
					max = TWO_D_ARRAY[b][r];
				if(min >TWO_D_ARRAY[b][r])	// Look for Min
					min = TWO_D_ARRAY[b][r];
			}
		
		}

		// Store Path Result
		RF_CH_Names[i] = FormatResult;		
		
		// Do a quick absolute value. We want the greatest distance from 0.0
		// to be the Max value.  We also want all data to be positive to reduce
		// the amount of work/time it takes to process the resulting data using
		// realtime analysis. (realtime has an absolute value function but adds
		// time to use. This includes string proccessing + math so I guess we
		// decided to handle it this way.)
		if(max < 0.0)
			max = max * -1.0;
		if(min < 0.0)
			min = min * -1.0;

		// Store Max and Min Values
		if(max > min)
		{
			RF_CH_Data[c]	= max;// Max Data
		//	RF_CH_Data[c+1] = min;// Min Data
		}
		else
		{
			RF_CH_Data[c]	= min;// Max Data
		//	RF_CH_Data[c+1] = max;// Min Data
		}
		
	}
	
	// Done
	// delete the 2-d array
	for(i = 0; i < column; ++i) 
		delete [] TWO_D_ARRAY[i];

	delete [] TWO_D_ARRAY;

}
/// @}
/// @}

/// @addtogroup Report
/// @{
/// @addtogroup PaceReport Format Pace report (XML).
/// @{

/** @brief <B>PaceReport:\n</B>
  This creates the pacing data XML data output.  It uses both the RV_Pace_Results &\n
  A_Pace_Results Arrays.  These arrays are initialized in main to -1 with a size of\n
  6 elements.  '-1' should be seen as no data.  This index result for each is as follows:\n\n

  Index[0] = PW_Max\n
  Index[1] = PW_Min\n   
  Index[2] = Ampl_Max\n 
  Index[3] = Ampl_Min\n
  Index[4] = Interval_Max\n  
  Index[5] = Interval_Min\n\n  

  Again, this is going to be done quickly not pretty.     

@param[in,out] error Error Status Information
*/	
TString PaceReport(ErrorCluster *error)
{
	//  This is a little different then normal.  This data is going to the results
	//  file.  If there was and error, the data may not be valid so I'm going to
	//  re-init the arrays to contain '-1.0'.  This is because I don't want to screw
	//  up the result file by having columns of different sizes.  This function must 
	//  return results so no errors can be created by this function.
	if (error->Status())
	{
		RV_Pace_Results.Initialize(-1.0, 6);	
		A_Pace_Results.Initialize(-1.0, 6);
	}

	TString DataOut = TString();			// Output Result
	TString RV = TString("RV");				// Channel Name
	TString Atrial = TString("Atrial");		// Channel Name
	TString FormatResult;

	// PW Results Defined
	TString PW_A = TString("<result_element>\n\t\t\t<name>");
	TString PW_B = TString("_PW_Max</name>\n\t\t\t<units>mSec</units>\n\t\t\t<type>SGL</type>\n\t\t\t<value>");
	TString PW_C = TString("</value>\n\t\t</result_element>\n<result_element>\n\t\t\t<name>");
	TString PW_D = TString("_PW_Min</name>\n\t\t\t<units>mSec</units>\n\t\t\t<type>SGL</type>\n\t\t\t<value>");
	TString PW_E = TString("</value>\n\t\t</result_element>\n");

	// Ampl Results Defined
	TString Ampl_A = TString("<result_element>\n\t\t\t<name>");
	TString Ampl_B = TString("_Ampl_Max</name>\n\t\t\t<units>Volts</units>\n\t\t\t<type>SGL</type>\n\t\t\t<value>");
	TString Ampl_C = TString("</value>\n\t\t</result_element>\n<result_element>\n\t\t\t<name>");
	TString Ampl_D = TString("_Ampl_Min</name>\n\t\t\t<units>Volts</units>\n\t\t\t<type>SGL</type>\n\t\t\t<value>");
	TString Ampl_E = TString("</value>\n\t\t</result_element>\n");

	// Interval Results Defined
	TString Interval_A = TString("<result_element>\n\t\t\t<name>");
	TString Interval_B = TString("_Interval_Max</name>\n\t\t\t<units>mSec</units>\n\t\t\t<type>SGL</type>\n\t\t\t<value>");
	TString Interval_C = TString("</value>\n\t\t</result_element>\n<result_element>\n\t\t\t<name>");
	TString Interval_D = TString("_Interval_Min</name>\n\t\t\t<units>mSec</units>\n\t\t\t<type>SGL</type>\n\t\t\t<value>");
	TString Interval_E = TString("</value>\n\t\t</result_element>\n");

	// Build the results (Start with RV)
	DataOut += PW_A;		// PW
	DataOut += RV;
	DataOut += PW_B;
	FormatResult.Format("%.6f", RV_Pace_Results[0]);
	DataOut += FormatResult;
	DataOut += PW_C;
	DataOut += RV;
	DataOut += PW_D;
	FormatResult.Format("%.6f", RV_Pace_Results[1]);
	DataOut += FormatResult;
	DataOut += PW_E;

	DataOut += Ampl_A;		// Amplitude
	DataOut += RV;
	DataOut += Ampl_B;
	FormatResult.Format("%.6f", RV_Pace_Results[2]);
	DataOut += FormatResult;
	DataOut += Ampl_C;
	DataOut += RV;
	DataOut += Ampl_D;
	FormatResult.Format("%.6f", RV_Pace_Results[3]);
	DataOut += FormatResult;
	DataOut += Ampl_E;

	DataOut += Interval_A;	// Interval
	DataOut += RV;
	DataOut += Interval_B;
	FormatResult.Format("%.6f", RV_Pace_Results[4]);
	DataOut += FormatResult;
	DataOut += Interval_C;
	DataOut += RV;
	DataOut += Interval_D;
	FormatResult.Format("%.6f", RV_Pace_Results[5]);
	DataOut += FormatResult;
	DataOut += Interval_E;
	
	// Build Atrial results 
	DataOut += PW_A;		// PW
	DataOut += Atrial;
	DataOut += PW_B;
	FormatResult.Format("%.6f", A_Pace_Results[0]);
	DataOut += FormatResult;
	DataOut += PW_C;
	DataOut += Atrial;
	DataOut += PW_D;
	FormatResult.Format("%.6f", A_Pace_Results[1]);
	DataOut += FormatResult;
	DataOut += PW_E;
	
	DataOut += Ampl_A;		// Amplitude
	DataOut += Atrial;
	DataOut += Ampl_B;
	FormatResult.Format("%.6f", A_Pace_Results[2]);
	DataOut += FormatResult;
	DataOut += Ampl_C;
	DataOut += Atrial;
	DataOut += Ampl_D;
	FormatResult.Format("%.6f", A_Pace_Results[3]);
	DataOut += FormatResult;
	DataOut += Ampl_E;
	
	DataOut += Interval_A;	// Interval
	DataOut += Atrial;
	DataOut += Interval_B;
	FormatResult.Format("%.6f", A_Pace_Results[4]);
	DataOut += FormatResult;
	DataOut += Interval_C;
	DataOut += Atrial;
	DataOut += Interval_D;
	FormatResult.Format("%.6f", A_Pace_Results[5]);
	DataOut += FormatResult;
	DataOut += Interval_E;

	return DataOut;

}
/// @}
/// @}

/// @addtogroup Report
/// @{
/// @addtogroup RF_Report Format RF report (XML).
/// @{

/** @brief <B>RF_Report:\n</B>
  This will be the rf rectification data report.  Two arrays are
  used to create this report.  RF_CH_Names contains all of the
  channel name combinations and RF_CH_Data contains the Max of the
  data collected.  This makes generating this report pretty
  easily done using a for loop.

@param[in,out] error Error Status Information
*/
TString RF_Report(ErrorCluster *error)
{
	
	// I shouldn't have gotten here if I didn't have any channel names.
	// If an error occured, I'm still going to create a report but I
	// will not use the channel data because I don't trust it.  I will
	// re-init the data array with a -1.  This should indicate no data
	// taken.
	if (error->Status())	
		RF_CH_Data.Initialize(-1.0, RF_CH_Names.GetSize()*2);

	TString DataOut = TString();	// Output Result
	TString FormatResult;
	int d = 0;
	
	// Define the sections of the result XML used
	TString Data_A = TString("<result_element>\n\t\t\t<name>");
	TString Data_B = TString("_Max</name>\n\t\t\t<units>Volts</units>\n\t\t\t<type>SGL</type>\n\t\t\t<value>");
	TString Data_C = TString("</value>\n\t\t</result_element>\n");

	// Build the result String
	for(int i = 0; i < RF_CH_Names.GetSize(); i++, d++)
	{
		DataOut += Data_A;
		DataOut += RF_CH_Names[i];
		DataOut += Data_B;
		FormatResult.Format("%.6f", RF_CH_Data[d]);	// Max Data
		DataOut += FormatResult;
		DataOut += Data_C;
	}

	return DataOut;

}
/// @}
/// @}

// FILTER INFORMATION STRUCTURE FOR FILTER ROUTINES

typedef struct {
	unsigned int length;	// size of filter
	float *history;			// pointer to history in filter
	float *coef;			// pointer to coefficients of filter
} FILTER;

#define FILTER_SECTIONS	 2	// 2 filter sections for 24 db/oct filter

typedef struct {
		double a0, a1, a2;       // numerator coefficients
		double b0, b1, b2;       // denominator coefficients
} BIQUAD;

void szxform(
	double *a0, double *a1, double *a2,	// numerator coefficients
	double *b0, double *b1, double *b2,	// denominator coefficients
	double fc,		// Filter cutoff frequency
	double fs,		// sampling rate
	double *k,		// overall gain factor
	float *coef);	// pointer to 4 iir coefficients

/// @addtogroup Filter
/// @{
/// @addtogroup iir_filter Perform IIR filtering sample by sample on floats.
/// @{

/** @brief <B>iir_filter:\n</B>
 iir_filter - Perform IIR filtering sample by sample on floats
 
 Implements cascaded direct form II second order sections.
 Requires FILTER structure for history and coefficients.
 The length in the filter structure specifies the number of sections.
 The size of the history array is 2*iir->length.
 The size of the coefficient array is 4*iir->length + 1 because
 the first coefficient is the overall scale factor for the filter.
 Returns one output sample for each input sample.  Allocates history
 array if not previously allocated.\n\n
 
 Returns float value giving the current output.\n

@param[in,out] error Error Status Information
@param[in] iir        pointer to FILTER structure
@param[in] input        new float input sample
*/					
float iir_filter(ErrorCluster *error, FILTER &iir, float input)
{
	if (error->Status())	// if error in, bail out
		return(0.00);

	unsigned int i = 0;
	float *hist1_ptr = NULL, *hist2_ptr = NULL, *coef_ptr = NULL;
	float output = 0, new_hist = 0, history1 = 0, history2 = 0;

	// allocate history array
	if(!iir.history)
		iir.history = (float *) calloc(2*iir.length,sizeof(float));

	if(!iir.history) 
	{
		error->SetError( -1, "iir_filter: Unable to allocate history");
		return(0.00);
	}

	coef_ptr = iir.coef;		// coefficient pointer
	hist1_ptr = iir.history;	// first history
	hist2_ptr = hist1_ptr + 1;	// next history

	// 1st number of coefficients array is overall input scale factor, or filter gain
	output = input * (*coef_ptr++);

	for (i = 0 ; i < iir.length; i++)
	{
		history1 = *hist1_ptr; // history values
		history2 = *hist2_ptr;
		
		output = output - history1 * (*coef_ptr++);
		new_hist = output - history2 * (*coef_ptr++); // poles

		output = new_hist + history1 * (*coef_ptr++);
		output = output + history2 * (*coef_ptr++);   // zeros

		*hist2_ptr++ = *hist1_ptr;
		*hist1_ptr++ = new_hist;

		hist1_ptr++;
		hist2_ptr++;
	}

	return(output);
}

/// @}
/// @}

/// @addtogroup Filter
/// @{
/// @addtogroup LowPass_Filter Perform Low Pass filter function.
/// @{

/** @brief <B>LowPass_Filter:\n</B>
Perform Low Pass filter function on data array handed to this function.  It will Update filter 
coefficients and create a 4th order filter (24 db/oct rolloff), consisting of two second order 
sections.  The filtered data will be returned.


@param[in,out] error Error Status Information
@param[in] fSample	Sample rate
@param[in] fCutOff	Cutoff frequency
@param[in] InData	Data Array being filtered
*/					
float32Array LowPass_Filter(ErrorCluster *error, float64 fSample, float64 fCutOff, float32Array InData)  
{
	if (error->Status())	// if error in, bail out
	  return float32Array();

	if(fCutOff > (fSample/2.00) - 1.00)
	{
		error->SetError(-1, "LP_Filter: Cutoff frequency > '%.02f'Hz not allowed.",(fSample/2.00) - 1.00);
		return float32Array();
	}

	// Setup filter s-domain coefficients
	BIQUAD ProtoCoef[FILTER_SECTIONS];

	// Section 1
	ProtoCoef[0].a0 = 1.0;
	ProtoCoef[0].a1 = 0;
	ProtoCoef[0].a2 = 0;
	ProtoCoef[0].b0 = 1.0;
	ProtoCoef[0].b1 = 0.765367;
	ProtoCoef[0].b2 = 1.0;

	// Section 2
	ProtoCoef[1].a0 = 1.0;
	ProtoCoef[1].a1 = 0;
	ProtoCoef[1].a2 = 0;
	ProtoCoef[1].b0 = 1.0;
	ProtoCoef[1].b1 = 1.847759;
	ProtoCoef[1].b2 = 1.0;

	FILTER   iir;
	iir.length = FILTER_SECTIONS;	// Number of filter sections
	iir.history = NULL;				// Initialize pointers
	iir.coef = NULL;

	// Allocate array of z-domain coefficients for each filter section
	// plus filter gain variable
	iir.coef = (float *) calloc(4 * iir.length + 1, sizeof(float));

	if (!iir.coef) // Shouldn't happen
	{
		error->SetError(-1, "LowPass Filter: Unable to allocate coef array.");
		return float32Array();
	}

	double k = 1.0;				// overall filter gain 
	float *coef = iir.coef + 1;	// Skip k, or gain 

	double Q = 1;				// Resonance 
	double fc = fCutOff;		// Filter cutoff (Hz)
	double fs = fSample;		// Sampling frequency (Hz) 

	// Compute z-domain coefficients for each biquad section
	// for new Cutoff Frequency and Resonance
	double   a0, a1, a2, b0, b1, b2;

	for( unsigned nInd = 0; nInd < iir.length; nInd++)
	{
		a0 = ProtoCoef[nInd].a0;
		a1 = ProtoCoef[nInd].a1;
		a2 = ProtoCoef[nInd].a2;

		b0 = ProtoCoef[nInd].b0;
		b1 = ProtoCoef[nInd].b1 / Q;	// Divide by resonance or Q
		b2 = ProtoCoef[nInd].b2;
		szxform(&a0, &a1, &a2, &b0, &b1, &b2, fc, fs, &k, coef);
		coef += 4;						// Point to next filter section 
	}

	// Update overall filter gain in coef array 
	iir.coef[0] = k;

	float32Array Conv_Data( InData.GetSize() );
	
	for(int z = 0; z < InData.GetSize(); z++)
		Conv_Data[z] = iir_filter(error, iir, InData[z]);

	free(iir.history);
	free(iir.coef);

	return Conv_Data;
}
/// @}
/// @}

/// @addtogroup Filter
/// @{
/// @addtogroup prewarp Perform bilinear transformation on s-domain coefficients.
/// @{

/** @brief <B>prewarp:\n</B>

      Perform bilinear transformation on s-domain coefficients\n
      of 2nd order biquad section.\n
      First design an analog filter and use s-domain coefficients\n
      as input to szxform() to convert them to z-domain.\n\n

      Here's the butterworth polynomials for 2nd, 4th and 6th order sections.\n
      When we construct a 24 db/oct filter, we take to 2nd order\n
      sections and compute the coefficients separately for each section.\n\n

      n       Polynomials\n
 --------------------------------------------------------------------\n
      2       s^2 + 1.4142s +1\n
      4       (s^2 + 0.765367s + 1) (s^2 + 1.847759s + 1)\n
      6       (s^2 + 0.5176387s + 1) (s^2 + 1.414214 + 1) (s^2 + 1.931852s + 1)\n\n

      Where n is a filter order.\n
      For n=4, or two second order sections, we have following equations for each\n
      2nd order stage:\n\n

      (1 / (s^2 + (1/Q) * 0.765367s + 1)) * (1 / (s^2 + (1/Q) * 1.847759s + 1))\n\n

      Where Q is filter quality factor in the range of\n
      1 to 1000. The overall filter Q is a product of all\n
      2nd order stages. For example, the 6th order filter\n
      (3 stages, or biquads) with individual Q of 2 will\n
      have filter Q = 2 * 2 * 2 = 8.\n\n

      The nominator part is just 1.\n
      The denominator coefficients for stage 1 of filter are:\n
      b2 = 1; b1 = 0.765367; b0 = 1;\n
      numerator is\n
      a2 = 0; a1 = 0; a0 = 1;\n\n

      The denominator coefficients for stage 1 of filter are:\n
      b2 = 1; b1 = 1.847759; b0 = 1;\n
      numerator is\n
      a2 = 0; a1 = 0; a0 = 1;\n\n

      These coefficients are used directly by the szxform()\n
      and bilinear() functions. For all stages the numerator\n
      is the same and the only thing that is different between\n
      different stages is 1st order coefficient. The rest of\n
      coefficients are the same for any stage and equal to 1.\n

*/					
void prewarp(double *a0, double *a1, double *a2, double fc, double fs);
void bilinear(
	double a0, double a1, double a2,	// numerator coefficients
	double b0, double b1, double b2,	// denominator coefficients
	double *k,							// overall gain factor
	double fs,							// sampling rate
	float *coef);						// pointer to 4 iir coefficients

void prewarp(
	double *a0, double *a1, double *a2,
	double fc, double fs)
{
	double wp, pi;

	pi = 4.0 * atan(1.0);
	wp = 2.0 * fs * tan(pi * fc / fs);

	*a2 = (*a2) / (wp * wp);
	*a1 = (*a1) / wp;
}

/// @}
/// @}

/// @addtogroup Filter
/// @{
/// @addtogroup bilinear Transform the numerator and denominator coefficients.
/// @{

/** @brief <B>bilinear:\n</B>

 Transform the numerator and denominator coefficients\n
 of s-domain biquad section into corresponding\n
 z-domain coefficients.\n\n

      Store the 4 IIR coefficients in array pointed by coef\n
      in following order:\n
             beta1, beta2    (denominator)\n
             alpha1, alpha2  (numerator)\n\n

 Arguments:\n
             a0-a2   - s-domain numerator coefficients\n\n
             b0-b2   - s-domain denominator coefficients\n\n
             k       - filter gain factor. initially set to 1\n
                       and modified by each biquad section in such\n
                       a way, as to make it the coefficient by\n
                       which to multiply the overall filter gain\n
                       in order to achieve a desired overall filter gain,\n
                       specified in initial value of k.\n\n
             fs      - sampling rate (Hz)\n\n
             coef    - array of z-domain coefficients to be filled in.\n\n

 Return:\n
             On return, set coef z-domain coefficients
*/
void bilinear(
	double a0, double a1, double a2,  // numerator coefficients
	double b0, double b1, double b2,  // denominator coefficients
	double *k,    // overall gain factor
	double fs,    // sampling rate
	float *coef   // pointer to 4 iir coefficients
)
{
	double ad, bd;

	// alpha (Numerator in s-domain)
	ad = 4.0 * a2 * fs * fs + 2.0 * a1 * fs + a0;
	// beta (Denominator in s-domain)
	bd = 4.0 * b2 * fs * fs + 2.0 * b1* fs + b0;

	// update gain constant for this section
	*k *= ad/bd;

	// Denominator
	*coef++ = (2.0 * b0 - 8.0 * b2 * fs * fs) / bd;				// beta1
	*coef++ = (4.0 * b2 * fs * fs - 2.0 * b1 * fs + b0) / bd;	// beta2

	// Nominator
	*coef++ = (2.0 * a0 - 8.0 * a2 * fs * fs) / ad;				// alpha1
	*coef = (4.0 * a2 * fs * fs - 2.0 * a1 * fs + a0) / ad;		// alpha2
}
/// @}
/// @}

/// @addtogroup Filter
/// @{
/// @addtogroup szxform Transform from s to z domain using bilinear transform with prewarp.
/// @{

/** @brief <B>szxform:\n</B>

  Transform from s to z domain using bilinear transform with prewarp.\n\n
 
  Arguments:\n
       For argument description look at bilinear()\n\n
 
       coef - pointer to array of floating point coefficients,\n
                      corresponding to output of bilinear transform\n
                      (z domain).\n\n
 
  Note: frequencies are in Hz.
*/
void szxform(
	double *a0, double *a1, double *a2,	// numerator coefficients
	double *b0, double *b1, double *b2,	// denominator coefficients
	double fc,		// Filter cutoff frequency
	double fs,		// sampling rate
	double *k,		// overall gain factor
	float *coef)	// pointer to 4 iir coefficients
{
	// Calculate a1 and a2 and overwrite the original values
	prewarp(a0, a1, a2, fc, fs);
	prewarp(b0, b1, b2, fc, fs);
	bilinear(*a0, *a1, *a2, *b0, *b1, *b2, k, fs, coef);
}
/// @}
/// @}

#pragma warning( pop )	// Restore Warning ('identifier' : has bad storage class)
