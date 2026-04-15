# Processing NHANES 2003-2006 semi-raw accelerometer data (ActiGraph 7164, single-axis, 60-second epoch counts) using the same GGIR pipeline as the 2011-2014 cycles for cross-cycle comparability.

## Background

The NHANES 2003 and 2005 cycles released semi-raw accelerometer data nearly two decades ago:
- [2003 Cycle (PAXRAW_C)](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/PAXRAW_C.htm)
- [2005 Cycle (PAXRAW_D)](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/PAXRAW_D.htm)

NCI has already released an [algorithm](https://epi.grants.cancer.gov/nhanes-pam/create.html) for processing this data, and there is a [great R package](https://github.com/vandomed/nhanesaccel) available as well. But my goal is different: I already processed the 2011-2014 cycles with GGIR, and I want results that are comparable across cycles using the same pipeline.

## Important Limitations

The 2003-2006 data were collected using the ActiGraph 7164, a single-axis accelerometer with 60-second epoch count data. Participants were instructed to remove the device while showering, swimming, or sleeping ([protocol](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/PAXRAW_C.htm)).

As a result:
1. Most sleep detection and step count algorithms (e.g., Verisense) do not work on this data.
2. Circadian rhythm estimates are unreliable because wear time is generally insufficient. L5 is most likely to coincide with the non-wear period (during sleep), producing RA values near 1 in ~90% of cases, which is not meaningful. Several published studies have reported RA, M10, and L5 from the 2003-2006 NHANES cycles. These results should be interpreted with caution, as the hip-worn, wake-only wear protocol means L5 almost always falls in the non-wear (sleep) period, artificially inflating RA toward 1 regardless of the participant's actual rest-activity pattern.
3. In roughly 10% of participants who did not follow the removal protocol, circadian rhythm parameters may be estimable.

## Step 0: Download

Download [PAXRAW_C](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/PAXRAW_C.zip) and [PAXRAW_D](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/PAXRAW_D.zip).

## Step 1: NCI Calibration

Run the NCI calibration using [create.pam_perminute.sas](https://epi.grants.cancer.gov/nhanes-pam/create.pam_perminute.sas). NCI appears to have manually reviewed all 2003-2004 data and performed QC (see lines 122-156 of the SAS code). I also recommend excluding records with PAXSTAT = 2.

This SAS file only covers 2003-2004. You can modify it for 2005-2006 by skipping the manual QC section. After this step you will have two SAS files: `pam_perminute_c.sas7bdat` and `pam_perminute_d.sas7bdat`, totaling about 8.43 GB.

## Step 2: Merge with Demographics

Merge age information from the DEMO files, since NCI uses age-specific cutoffs (see [create.pam_perday.sas](https://epi.grants.cancer.gov/nhanes-pam/create.pam_perday.sas), lines 206-250).

```r
pam_perminute_C <- read_sas("pam_perminute_c.sas7bdat")
pam_perminute_D <- read_sas("pam_perminute_d.sas7bdat")
demo_C <- read_xpt("DEMO_C.XPT")
demo_D <- read_xpt("DEMO_D.XPT")
demo <- bind_rows(demo_C, demo_D) %>% dplyr::select(SEQN, RIDAGEYR)
pam_perminute <- bind_rows(pam_perminute_C, pam_perminute_D %>% dplyr::select(-PAXSTEP))
pam_perminute <- inner_join(pam_perminute, demo, by = "SEQN")
```

## Step 3: Generate ActiGraph-Style CSV Files

Format the data into ActiGraph GT3X-style CSV files using VM = count, matching the format used in the 2011-2014 pipeline. This ensures comparability and avoids writing a separate reader.

The code below processes adults (age >= 18). Age-specific cutoffs require separate processing for other age groups.

<details>
<summary>CSV generation code</summary>

```r
pam_perminuteWorking <- filter(pam_perminute, RIDAGEYR >= 18) %>%
  dplyr::select(-RIDAGEYR)
allSEQN <- unique(pam_perminuteWorking$SEQN)

# Map PAXDAY (1=Sunday...7=Saturday) to a date
first_week <- as.Date("2000-01-01") + 0:6
day_codes  <- as.POSIXlt(first_week)$wday + 1
date_map   <- setNames(first_week, day_codes)

for (seq in allSEQN) {
  temp <- filter(pam_perminuteWorking, SEQN == seq)

  out <- temp %>%
    mutate(
      Date = date_map[as.character(PAXDAY)],
      TimeStamp = as.POSIXct(
        paste0(Date, " ", sprintf("%02d:%02d:00", PAXHOUR, PAXMINUT)),
        tz = "UTC"
      ),
      TimeStamp = format(TimeStamp, "%Y-%m-%dT%H:%M:%SZ"),
      vm    = PAXINTEN,
      axis1 = vm, axis2 = 0, axis3 = 0  # dummy values; GGIR uses VM only
    ) %>%
    select(TimeStamp, axis1, axis2, axis3, vm)

  hdr <- c(
    "------------ Data Table File Created By Actigraph Link ActiLife v6.11.9 date format dd/MM/yyyy Filter Normal -----------,,,,,",
    "Serial Number: TAS1D48140206,,,,,",
    paste0("Start Time ",
      format(as.POSIXct(paste0(date_map[as.character(temp$PAXDAY[1])], " ",
        sprintf("%02d:%02d:00", temp$PAXHOUR[1], temp$PAXMINUT[1])), tz = "UTC"), "%H:%M:%S"), ",,,,,"),
    paste0("Start Date ",
      format(as.Date(date_map[as.character(temp$PAXDAY[1])]), "%m-%d-%Y"), ",,,,,"),
    "Epoch Period (hh:mm:ss) 00:00:60,,,,,",
    "Download Time 15:29:44,,,,,",
    "Download Date 09-19-2017,,,,,",
    "Current Memory Address: 0,,,,,",
    "Current Battery Voltage: 3.85     Mode = 13,,,,,",
    "--------------------------------------------------,,,,,"
  )

  file_path <- file.path("/NHANESGGIR/NHANES03040506Adult/", paste0(seq, ".csv"))
  writeLines(hdr, file_path)
  write.table(out, file_path, sep = ",", row.names = FALSE,
              col.names = TRUE, quote = FALSE, append = TRUE)
}
```

</details>

## Step 4: Process with GGIR

Key design decisions for using GGIR on 2003-2006 data:

- **No calibration or imputation**: not feasible with single-axis count data (`do.cal = FALSE`, `do.imp = FALSE`)
- **Window sizes**: `c(60, 3600, 3600)` to match the NCI algorithm (60s native epoch, 1-hour non-wear and assessment windows)
- **Non-wear as sleep proxy**: uses GGIR's NotWorn algorithm for both `HASPT.algo` and `HASIB.algo`, since actual sleep detection is not possible with this data
- **10-hour minimum wear time** per valid day, matching NCI
- **NCI activity thresholds**: sedentary < 100 cpm, MVPA >= 2020 cpm (age-specific, adjust as needed)
- **NeishabouriCount columns are not usable**: they contain dummy values from the CSV formatting step

<details>
<summary>GGIR configuration</summary>

```r
GGIR(
  idloc                    = 6,
  datadir                  = "/NHANESGGIR/NHANES03040506Adult",
  outputdir                = "/NHANESGGIR/ResultNHANES03040506Adult",
  do.cal                   = FALSE,
  dataFormat               = "actigraph_csv",
  extEpochData_timeformat  = "%m-%d-%Y %H:%M:%S",
  windowsizes              = c(60, 3600, 3600),
  do.neishabouricounts     = FALSE,
  acc.metric               = "NeishabouriCount_vm",
  studyname                = "NHANES03040506",
  HASPT.algo               = "NotWorn",
  HASIB.algo               = "NotWorn",
  do.imp                   = FALSE,
  nonwear_range_threshold  = 100,
  HASPT.ignore.invalid     = NA,
  ignorenonwear            = FALSE,
  threshold.lig            = 100,
  threshold.mod            = 2020,
  threshold.vig            = 5999,
  includedaycrit           = 10,
  includenightcrit         = 10
)
```

</details>

## Validation

I compared MVPA and sedentary time against the NCI SAS algorithm for participants aged >= 18 (raw correlation, no survey weights):

| Metric | Correlation with NCI |
|---|---|
| MVPA | 97.12% |
| Sedentary time | 85.43% |

MVPA correlation is high because at the 60-second epoch level, bout detection is unnecessary and the count threshold is straightforward. The lower sedentary correlation likely reflects differences in non-wear handling between GGIR and the NCI algorithm.

## Contact

If you have questions or would like access to the processed data (too large for GitHub), feel free to reach out on [LinkedIn](https://www.linkedin.com/in/xiaoyu-zong-0a733ba0/)
