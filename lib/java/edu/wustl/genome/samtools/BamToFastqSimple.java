package edu.wustl.genome.samtools;

import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import net.sf.picard.io.IoUtil;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMUtils;
import net.sf.samtools.util.SequenceUtil;
import net.sf.samtools.util.StringUtil;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;



public class BamToFastqSimple {

    private String bamFileName;
    private File INPUT;
    private final Map<String,SAMRecord> firstSeenMates;
    private final SAMFileReader reader;
    private final SAMRecordIterator recordIterator;

    public static void main(String[] args) {
        BamToFastqSimple b = new BamToFastqSimple(args[0]); 
        try {
            b.doWork();
        } catch (Exception e) {
            System.out.println("oh nos.  something went wrong: " + e);
        }
    }

    public BamToFastqSimple(String bamFileName) {
       this.bamFileName = bamFileName;
       INPUT = new File(bamFileName);  
       firstSeenMates = new HashMap<String,SAMRecord>();
       reader = new SAMFileReader(IoUtil.openFileForReading(INPUT));
       recordIterator = reader.iterator();
    }

    public void finished() throws Exception {
        reader.close();
        if (firstSeenMates.size() > 0) {
            throw new Exception("Found "+firstSeenMates.size()+" unpaired mates");
        }
    }

    public boolean hasRemainingData() {
        return recordIterator.hasNext();
    }

    public void doWork() throws Exception {
      while (hasRemainingData())
         System.out.println(getNextFastqPair()); 
    }

    public String getNextFastqPair() throws Exception {
        while (recordIterator.hasNext()) {
            SAMRecord currentRecord = recordIterator.next();
            // Skip non-PF reads as necessary
            if (currentRecord.getReadFailsVendorQualityCheckFlag()) continue;

            final String currentReadName = currentRecord.getReadName() ;
            final SAMRecord firstRecord = firstSeenMates.get(currentReadName);
            if (firstRecord == null) {
                firstSeenMates.put(currentReadName, currentRecord) ;
            }
            else {
                assertPairedMates(firstRecord, currentRecord);
                String ret;
                if (currentRecord.getFirstOfPairFlag()) {
                     ret = samRecordToFastq(currentRecord, 1) + "\n" + 
                            samRecordToFastq(firstRecord, 2);
                }
                else {
                     ret = samRecordToFastq(firstRecord, 1) + "\n" + 
                            samRecordToFastq(currentRecord, 2);
                }
                firstSeenMates.remove(currentReadName);
                return ret;
            }
        }
        return null;
    }

    private void assertPairedMates(final SAMRecord record1, final SAMRecord record2) throws Exception {
        if (! (record1.getFirstOfPairFlag() && record2.getSecondOfPairFlag() ||
               record2.getFirstOfPairFlag() && record1.getSecondOfPairFlag() ) ) {
            throw new Exception("Illegal mate state: " + record1.getReadName());
        }
    }

    private String samRecordToFastq(final SAMRecord read, final Integer mateNumber) {
        final String seqHeader = mateNumber==null ? read.getReadName() : read.getReadName() + "/"+ mateNumber;
        String readString = read.getReadString();
        String baseQualities = read.getBaseQualityString();

        // by default reverse complement things here until we decide otherwise
        if ( read.getReadNegativeStrandFlag() ) {
            readString = SequenceUtil.reverseComplement(readString);
            baseQualities = StringUtil.reverseString(baseQualities);
        }

        return "@" + seqHeader + "\n" + readString + "\n+\n" + baseQualities;
    }
}
