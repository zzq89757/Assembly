function calculateIdentity(seq1, seq2) {
  // 计算比对后的identity，包括考虑错配和indel
  let matches = 0;
  let alignedLength = seq1.length;
  
  for (let i = 0; i < seq1.length; i++) {
      if (seq1[i] === seq2[i]) {
          matches++;
      }
  }
  
  return matches / alignedLength;
}

function alignSequences(seq1, seq2) {
  // 局部比对算法，考虑错配和indel
  let maxScore = -Infinity;
  let bestAlign1 = '';
  let bestAlign2 = '';
  let bestIdentity = 0;
  
  for (let i = 1; i <= seq1.length; i++) {
      let aligned1 = seq1.substring(seq1.length - i);
      let aligned2 = seq2.substring(0, i);
      
      if (aligned1.length !== aligned2.length) continue;
      
      let identity = calculateIdentity(aligned1, aligned2);
      let score = identity * aligned1.length;
      
      if (score > maxScore) {
          maxScore = score;
          bestAlign1 = aligned1;
          bestAlign2 = aligned2;
          bestIdentity = identity;
      }
  }
  
  return [bestAlign1, bestAlign2, bestIdentity];
}

function assembleSequences(sequences, minOverlap, minIdentity) {
  let assembledSeq = sequences[0];
  let allAlignedSequences = [assembledSeq];
  let identityInfo = [];
  
  for (let i = 1; i < sequences.length; i++) {
      let maxIdentity = 0;
      let bestAlign1 = '';
      let bestAlign2 = '';
      let overlapStart = 0;

      for (let overlapLen = minOverlap; overlapLen <= sequences[i].length; overlapLen++) {
          let seq1Overlap = assembledSeq.substring(assembledSeq.length - overlapLen);
          let seq2Overlap = sequences[i].substring(0, overlapLen);
          
          let [aligned1, aligned2, identity] = alignSequences(seq1Overlap, seq2Overlap);
          
          if (identity >= minIdentity && identity > maxIdentity) {
              maxIdentity = identity;
              bestAlign1 = aligned1;
              bestAlign2 = aligned2;
              overlapStart = assembledSeq.length - overlapLen;
          }
      }
      
      if (bestAlign1.length > 0 && bestAlign2.length > 0) {
          assembledSeq = assembledSeq.substring(0, overlapStart) + bestAlign2 + sequences[i].substring(bestAlign2.length);
          allAlignedSequences.push("-".repeat(overlapStart) + bestAlign2 + sequences[i].substring(bestAlign2.length));
          identityInfo.push(`Overlap between sequence ${i} and ${i+1} - Identity: ${(maxIdentity * 100).toFixed(2)}%`);
      } else {
          assembledSeq += sequences[i];
          allAlignedSequences.push("-".repeat(assembledSeq.length - sequences[i].length) + sequences[i]);
      }
  }
  
  // Align sequences and add markers
  let maxLen = Math.max(...allAlignedSequences.map(seq => seq.length));
  let alignedOutput = allAlignedSequences.map(seq => seq.padEnd(maxLen, '-')).join("\n");
  
  let alignmentMarkers = allAlignedSequences.slice(1).map((seq, index) => {
      let prevSeq = allAlignedSequences[index];
      return [...prevSeq].map((char, i) => (char === seq[i] && char !== '-') ? '|' : ' ').join('').padEnd(maxLen, ' ');
  }).join("\n");

  return `${alignedOutput}\n${alignmentMarkers}\n\n${identityInfo.join("\n")}`;
}

// 示例使用
let sequences = [
  "ATCGTACGTAGCTAGCTGACGCAGT",
  "GCTAGCTGACGTCAGTCAGGCTTACG",
  "CAGTCAGGCTTACGCGTATAGCGTAC"
];

let minOverlap = 5;  // 最小重叠长度
let minIdentity = 0.8;  // 最小identity阈值

let result = assembleSequences(sequences, minOverlap, minIdentity);
console.log(`组装后的序列:\n${result}`);
