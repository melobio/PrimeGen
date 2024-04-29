/*
 Navicat Premium Data Transfer

 Source Server         : 120.24.193.47
 Source Server Type    : MySQL
 Source Server Version : 50736
 Source Host           : 120.24.193.47:3306
 Source Schema         : mgi_primer

 Target Server Type    : MySQL
 Target Server Version : 50736
 File Encoding         : 65001

 Date: 22/02/2024 14:14:34
*/

SET NAMES utf8mb4;
SET FOREIGN_KEY_CHECKS = 0;

-- ----------------------------
-- Table structure for dimer_info
-- ----------------------------
CREATE TABLE IF NOT EXISTS `dimer_info` (
  `id` bigint(20) unsigned NOT NULL AUTO_INCREMENT COMMENT 'primary id',
  `experiment_id` varchar(50) NOT NULL COMMENT 'experiment name',
  `experiment_order` varchar(50) DEFAULT NULL COMMENT 'Sample name (duplicate group)',
  `dimer` varchar(256) DEFAULT NULL COMMENT 'dimer sequence',
  `r1_primer` varchar(50) DEFAULT NULL COMMENT 'Primer 1 forming the dimer',
  `r2_primer` varchar(50) DEFAULT NULL COMMENT 'Primer 2 forming the dimer',
  PRIMARY KEY (`id`),
  KEY `dimer_info_experiment_id_IDX` (`experiment_id`) USING BTREE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COMMENT='The dimer sequence file of multiplex PCR provides the specific sequence of the dimer and the primer sequence that forms the dimer' ROW_FORMAT = Dynamic;
-- ----------------------------
-- Table structure for dt_wq_primer
-- ----------------------------
CREATE TABLE IF NOT EXISTS `dt_wq_primer`  (
  `id` bigint(20) UNSIGNED NOT NULL AUTO_INCREMENT COMMENT 'primary id',
  `experiment_id` varchar(50) CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci NOT NULL COMMENT 'experiment name',
  `name_id` varchar(50) CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci NULL DEFAULT NULL COMMENT 'Amplicon name',
  `chr_name1` varchar(50) CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci NULL DEFAULT NULL COMMENT 'Name of the gene sequence where primer1 of the amplicon is located',
  `primer1` varchar(256) CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci NULL DEFAULT NULL COMMENT 'Sequence of primer1 of the amplicon',
  `start1` int(11) NULL DEFAULT NULL COMMENT 'Starting position of primer1 of the amplicon',
  `end1` int(11) NULL DEFAULT NULL COMMENT 'Termination position of primer1 of the amplicon',
  `strand1` char(1) CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci NULL DEFAULT NULL COMMENT 'The chain on which the amplicon primer 1 is located',
  `chr_name2` varchar(50) CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci NULL DEFAULT NULL COMMENT 'Name of the gene sequence where primer 2 of the amplicon is located',
  `primer2` varchar(256) CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci NULL DEFAULT NULL COMMENT 'Sequence of primer2 of the amplicon',
  `start2` int(11) NULL DEFAULT NULL COMMENT 'Starting position of primer2 of the amplicon',
  `end2` int(11) NULL DEFAULT NULL COMMENT 'Termination position of primer2 of the amplicon',
  `strand2` char(1) CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci NULL DEFAULT NULL COMMENT 'The chain on which the amplicon primer2 is located',
  `create_time` datetime(0) NULL DEFAULT CURRENT_TIMESTAMP(0) COMMENT 'create time',
  PRIMARY KEY (`id`) USING BTREE,
  INDEX `dt_wq_primer_experiment_id_IDX`(`experiment_id`) USING BTREE
) ENGINE = InnoDB CHARACTER SET = utf8mb4 COLLATE = utf8mb4_general_ci COMMENT = 'Primer list file for multiplex PCR, providing sample amplicon names, primer sequences, primer start and terminate positions, sequence of the genome in which the primer is located, and the positive and negative chains of the primer in the genome' ROW_FORMAT = Dynamic;

-- ----------------------------
-- Table structure for multi_pcr_depth
-- ----------------------------
CREATE TABLE IF NOT EXISTS `multi_pcr_depth`  (
  `id` bigint(20) UNSIGNED NOT NULL AUTO_INCREMENT COMMENT 'primary id',
  `experiment_id` varchar(50) CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci NOT NULL COMMENT 'Experiment name',
  `experiment_order` varchar(50) CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci NULL DEFAULT NULL COMMENT 'Sample name (duplicate group)',
  `name` varchar(50) CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci NULL DEFAULT NULL COMMENT 'Amplicon name',
  `depth` decimal(10, 2) NULL DEFAULT NULL COMMENT 'average depth of coverage for the corresponding region in the sample',
  PRIMARY KEY (`id`) USING BTREE,
  INDEX `multiPcr_depth_experiment_id_IDX`(`experiment_id`) USING BTREE
) ENGINE = InnoDB CHARACTER SET = utf8mb4 COLLATE = utf8mb4_general_ci COMMENT = 'Sequencing depth analysis file for multiplex PCR, providing duplicate sample name, amplicon name and sequencing depth for analyzing the sequencing depth of amplicons under the sample' ROW_FORMAT = Dynamic;

-- ----------------------------
-- Table structure for multi_pcr_efficiency
-- ----------------------------
CREATE TABLE IF NOT EXISTS `multi_pcr_efficiency`  (
  `id` bigint(20) UNSIGNED NOT NULL AUTO_INCREMENT COMMENT 'primary id',
  `experiment_id` varchar(50) CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci NOT NULL COMMENT 'Experiment name',
  `experiment_order` varchar(50) CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci NULL DEFAULT NULL COMMENT 'Sample name (duplicate group)',
  `amp_name` varchar(50) CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci NULL DEFAULT NULL COMMENT 'Amplicon name',
  `f_bias` int(11) NULL DEFAULT NULL COMMENT 'Number of non-specific reads for the forward primer',
  `r_bias` int(11) NULL DEFAULT NULL COMMENT 'Number of non-specific reads for the reverse primer',
  `primer_efficiency` decimal(5, 2) NULL DEFAULT NULL COMMENT 'Primer amplification efficiency',
  PRIMARY KEY (`id`) USING BTREE,
  INDEX `multiPcr_efficiency_experiment_id_IDX`(`experiment_id`) USING BTREE
) ENGINE = InnoDB CHARACTER SET = utf8mb4 COLLATE = utf8mb4_general_ci COMMENT = 'Primer efficiency analysis file for multiplex PCR, providing amplicon names, and their corresponding number of forward and reverse primer non-specific reads and primer efficiency' ROW_FORMAT = Dynamic;

-- ----------------------------
-- Table structure for multi_pcr_summary
-- ----------------------------
CREATE TABLE IF NOT EXISTS `multi_pcr_summary`  (
  `id` bigint(20) UNSIGNED NOT NULL AUTO_INCREMENT COMMENT 'primary id',
  `experiment_id` varchar(50) CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci NOT NULL COMMENT 'Experiment name',
  `sample_name` varchar(50) CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci NULL DEFAULT NULL COMMENT 'Sample name (duplicate group)',
  `dimer_rate` decimal(5, 2) NULL DEFAULT NULL COMMENT 'The proportion of adapter-containing reads in the total reads, i.e., the dimer rate',
  `uniformity` decimal(5, 2) NULL DEFAULT NULL COMMENT 'Amplification homogeneity of all amplicons in the sample (replicate group)',
  `map_target_rate` decimal(5, 2) NULL DEFAULT NULL COMMENT 'The proportion of reads that specifically map to the target regions of the reference genome sequence, i.e., map target rate ',
  `cov_100x` decimal(5, 2) NULL DEFAULT NULL COMMENT 'Percentage of sequencing target regions sequenced more than or equal to 100 times, i.e., 100x depth ratio',
  PRIMARY KEY (`id`) USING BTREE,
  INDEX `summary_experiment_id_IDX`(`experiment_id`) USING BTREE
) ENGINE = InnoDB CHARACTER SET = utf8mb4 COLLATE = utf8mb4_general_ci COMMENT = 'Summary file for multiplex PCR analysis, providing sample name, sample dimer rate, sample amplification homogeneity, amplicon specific match rate for the sample, and sequencing depth' ROW_FORMAT = Dynamic;

-- ----------------------------
-- Table structure for multi_pcr_unexpected
-- ----------------------------
CREATE TABLE IF NOT EXISTS `multi_pcr_unexpected`  (
  `id` bigint(20) UNSIGNED NOT NULL AUTO_INCREMENT COMMENT 'primary id',
  `experiment_id` varchar(50) CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci NOT NULL COMMENT 'Experiment name',
  `experiment_order` varchar(50) CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci NULL DEFAULT NULL COMMENT 'Sample name (duplicate group)',
  `fprimer_rprimer` varchar(50) CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci NULL DEFAULT NULL COMMENT 'Forward Primer:Reverse Primer',
  `count` int(11) NULL DEFAULT NULL COMMENT 'Number of non-specific amplified reads',
  PRIMARY KEY (`id`) USING BTREE,
  INDEX `multiPcr_unexpected_experiment_id_IDX`(`experiment_id`) USING BTREE
) ENGINE = InnoDB CHARACTER SET = utf8mb4 COLLATE = utf8mb4_general_ci COMMENT = 'An analysis file for non-specific fragments of multiplex PCR, providing the sample name, forward and reverse primers, and the number of non-specific amplified reads, can be used to analyze non-specific amplified fragments.' ROW_FORMAT = Dynamic;

SET FOREIGN_KEY_CHECKS = 1;
