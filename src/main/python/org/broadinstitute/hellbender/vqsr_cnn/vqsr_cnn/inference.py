# Imports
import os
import sys
import h5py
import time
import numpy as np
from collections import Counter, defaultdict

# Keras Imports
import keras.backend as K
from keras.models import load_model

dna_symbols = {'A':0, 'C':1, 'G':2, 'T':3}

# Base calling ambiguities, See https://www.bioinformatics.org/sms/iupac.html
ambiguity_codes = {'K':[0,0,0.5,0.5], 'M':[0.5,0.5,0,0], 'R':[0.5,0,0,0.5], 'Y':[0,0.5,0.5,0], 'S':[0,0.5,0,0.5], 'W':[0.5,0,0.5,0], 
				  'B':[0,0.333,0.333,0.334], 'V':[0.333,0.333,0,0.334],'H':[0.333,0.333,0.334,0],'D':[0.333,0,0.333,0.334],
				  'X':[0.25,0.25,0.25,0.25], 'N':[0.25,0.25,0.25,0.25]}


annotations = ['MQ', 'DP', 'SOR', 'FS', 'QD', 'MQRankSum', 'ReadPosRankSum']

snp_indel_labels = {'NOT_SNP':0, 'NOT_INDEL':1, 'SNP':2, 'INDEL':3}
eps = 1e-7


def score_and_write_batch(model, file_out, fifo, batch_size, python_batch_size):
	'''Score and write a batch of variants with a 1D CNN.

	This function is tightly coupled with the NeuralNetInference
	It requires data written to the fifo in the order given by transferToPythonViaFifo

	Arguments
		model: a keras model
		file_out: The VCF file where variants scores are written
		fifo: The fifo opened by GATK Streaming executor
		batch_size: The total number of variants available in the fifo
		python_batch_size: the number of variants to process in each inference

	'''
	annotation_batch = []
	variant_data = []
	dna_batch = []
	score_key = ""
	genotypes = []
	is_snps = []

	for i in range(batch_size):
		fifo_line = fifo.readline()
		fifo_data = fifo_line.split('|')

		dna_batch.append(reference_string_to_tensor(fifo_data[0]))
		annotation_batch.append(annotation_string_to_tensor(fifo_data[1]))
		variant_data.append(fifo_data[2])
		score_key = fifo_data[3]
		genotypes.append(fifo_data[4])
		is_snps.append(int(fifo_data[5]))

	predictions = model.predict([np.array(dna_batch), np.array(annotation_batch)], batch_size=python_batch_size)
	snp_scores = predictions_to_snp_scores(predictions)
	indel_scores = predictions_to_indel_scores(predictions)

	for i in range(batch_size):
		if is_snps[i]:
			file_out.write(variant_data[i]+score_key+'={0:.3f}'.format(snp_scores[i])+genotypes[i]+'\n')
		else:
			file_out.write(variant_data[i]+score_key+'={0:.3f}'.format(indel_scores[i])+genotypes[i]+'\n')


def write_snp_score(model, file_out, fifo_line):
	fifo_data = fifo_line.split('|')
	score = snp_score_from_reference_annotations(model, fifo_data[0],  fifo_data[1])
	file_out.write(fifo_data[2] + fifo_data[3] + '=' + str(score) + fifo_data[4])


def write_indel_score(model, file_out, fifo_line):
	fifo_data = fifo_line.split('|')
	score = indel_score_from_reference_annotations(model, fifo_data[0],  fifo_data[1])
	file_out.write(fifo_data[2] + fifo_data[3] + '=' + str(score) + fifo_data[4])


def snp_score_from_reference_annotations(model, reference, annotation_string):
	dna_tensor = reference_string_to_tensor(reference)
	annotation_tensor = annotation_string_to_tensor(annotation_string)	
	
	dna_tensor = np.expand_dims(dna_tensor, axis=0)
	annotation_tensor = np.expand_dims(annotation_tensor, axis=0)

	predictions = model.predict([dna_tensor, annotation_tensor])[0]
	return np.log(eps + predictions[snp_indel_labels['SNP']] / (predictions[snp_indel_labels['NOT_SNP']] + eps))


def indel_score_from_reference_annotations(model, reference, annotation_string):
	dna_tensor = reference_string_to_tensor(reference)
	annotation_tensor = annotation_string_to_tensor(annotation_string)	
	
	dna_tensor = np.expand_dims(dna_tensor, axis=0)
	annotation_tensor = np.expand_dims(annotation_tensor, axis=0)

	predictions = model.predict([dna_tensor, annotation_tensor])[0]
	return np.log(eps + predictions[snp_indel_labels['INDEL']] / (predictions[snp_indel_labels['NOT_INDEL']] + eps))


def reference_string_to_tensor(reference):
	dna_data = np.zeros( (len(reference), len(dna_symbols)) )
	for i,b in enumerate(reference):
		if b in dna_symbols:
			dna_data[i, dna_symbols[b]] = 1.0
		elif b in ambiguity_codes:
			dna_data[i] = ambiguity_codes[b]
		else:
			print('Error! Unknown code:', b)
			continue
	return dna_data


def annotation_string_to_tensor(annotation_string):
	name_val_pairs = annotation_string.split(';')
	name_val_arrays = [p.split('=') for p in name_val_pairs]
	annotation_map = {str(p[0]).strip() : p[1] for p in name_val_arrays if len(p) > 1}
	annotation_data = np.zeros(( len(annotations), ))
	
	for i,a in enumerate(annotations):
		if a in annotation_map:
			annotation_data[i] = annotation_map[a]
	
	return annotation_data


def apply_model_to_batch(args, model, dna_batch, annotation_batch, variant_batch, vcf_writer, stats):
	predictions = model.predict([dna_batch, annotation_batch], batch_size=args.batch_size)
	snp_dict = predictions_to_snp_scores(args, predictions, positions)
	indel_dict = predictions_to_indel_scores(args, predictions, positions)

	# loop over the batch of variants and write them out with a score
	for v_out in variant_batch:
		position = v_out.contig + '_' + str(v_out.pos)

		if len(v_out.ref) == 1 and len(v_out.alleles[1][0]) == 1: # SNP means ref and alt both are length 1
			v_out.info['CNN_SCORE'] = float(snp_dict[position])
		else:
			v_out.info['CNN_SCORE'] = float(indel_dict[position])

		vcf_writer.write(v_out)
		stats['variants_written'] += 1


def predictions_to_snp_scores(predictions):
	snp = predictions[:, snp_indel_labels['SNP']]
	not_snp = predictions[:, snp_indel_labels['NOT_SNP']]
	return np.log(eps + snp / (not_snp + eps))


def predictions_to_indel_scores(predictions):
	indel = predictions[:, snp_indel_labels['INDEL']]
	not_indel = predictions[:, snp_indel_labels['NOT_INDEL']]
	return np.log(eps + indel / (not_indel + eps))


def predictions_to_snp_indel_scores(predictions):
	snp_dict = predictions_to_snp_scores(predictions)
	indel_dict = predictions_to_indel_scores(predictions)
	return snp_dict, indel_dict


def get_tensor_channel_map_1d_dna():
	'''1D Reference tensor with 4 channel DNA encoding.'''
	tensor_map = {}
	for k in inputs.keys():
		tensor_map[k] = inputs[k]
	
	return tensor_map


def get_tensor_channel_map_1d():
	'''1D Reference tensor with 4 channel DNA encoding'''
	tensor_map = {}
	for k in inputs.keys():
		tensor_map[k] = inputs[k]
	
	return tensor_map


def interval_file_to_dict(interval_file, shift1=0, skip=['@']):
	''' Create a dict to store intervals from a interval list file.

	Arguments:
		interval_file: the file to load either a bed file -> shift1 should be 1
			or a picard style interval_list file -> shift1 should be 0
		shift1: Shift the intervals 1 position over to align with 1-indexed VCFs
		skip: Comment character to ignore
	Returns:
		intervals: dict where keys in the dict are contig ids
			values are a tuple of arrays the first array 
			in the tuple contains the start positions
			the second array contains the end positions.
	'''
	intervals = {}

	with open(interval_file)as f:
		for line in f:
			if line[0] in skip:
				continue

			parts = line.split()
			contig = parts[0]
			lower = int(parts[1])+shift1
			upper = int(parts[2])+shift1

			if contig not in intervals:
				intervals[contig] = ([], [])

			intervals[contig][0].append(lower)
			intervals[contig][1].append(upper)

	for k in intervals.keys():
		intervals[k] = (np.array(intervals[k][0]), np.array(intervals[k][1]))		

	return intervals



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~ Metrics ~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def precision(y_true, y_pred):
	'''Calculates the precision, a metric for multi-label classification of
	how many selected items are relevant.
	'''
	true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
	predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
	precision = true_positives / (predicted_positives + K.epsilon())
	return precision


def recall(y_true, y_pred):
	'''Calculates the recall, a metric for multi-label classification of
	how many relevant items are selected.
	'''
	true_positives = K.sum(K.round(K.clip(y_true*y_pred, 0, 1)))
	possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
	recall = true_positives / (possible_positives + K.epsilon())
	return recall


def per_class_recall(labels):
	recall_fxns = []

	for label_key in labels:
		label_idx = labels[label_key]
		string_fxn = 'def '+ label_key + '_recall(y_true, y_pred):\n'
		string_fxn += '\ttrue_positives = K.sum(K.round(K.clip(y_true*y_pred, 0, 1)), axis=0)\n'
		string_fxn += '\tpossible_positives = K.sum(K.round(K.clip(y_true, 0, 1)), axis=0)\n'
		string_fxn += '\treturn true_positives['+str(label_idx)+'] / (possible_positives['+str(label_idx)+'] + K.epsilon())\n'

		exec(string_fxn)
		recall_fxn = eval(label_key + '_recall')
		recall_fxns.append(recall_fxn)

	return recall_fxns


def per_class_precision(labels):
	precision_fxns = []

	for label_key in labels:
		label_idx = labels[label_key]
		string_fxn = 'def '+ label_key + '_precision(y_true, y_pred):\n'
		string_fxn += '\ttrue_positives = K.sum(K.round(K.clip(y_true*y_pred, 0, 1)), axis=0)\n'
		string_fxn += '\tpredicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)), axis=0)\n'
		string_fxn += '\treturn true_positives['+str(label_idx)+'] / (predicted_positives['+str(label_idx)+'] + K.epsilon())\n'

		exec(string_fxn)
		precision_fxn = eval(label_key + '_precision')
		precision_fxns.append(precision_fxn)

	return precision_fxns


def get_metric_dict(labels=snp_indel_labels):
	metrics = {'precision':precision, 'recall':recall}
	precision_fxns = per_class_precision(labels)
	recall_fxns = per_class_recall(labels)
	for i,label_key in enumerate(labels.keys()):
		metrics[label_key+'_precision'] = precision_fxns[i]
		metrics[label_key+'_recall'] = recall_fxns[i]

	return metrics