export interface JetsonResult {
  event_name: 'bad_tip_change';
  pcr: number;
  tip: number;
  bad_tip: number;
  plot: Array<number>;
  data: string;
}
