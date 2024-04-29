import { createHash } from 'crypto'

export function calculateSHA256Hash(message: string) {
  const sha256Hash = createHash('sha256');
  sha256Hash.update(message);
  return sha256Hash.digest('hex');
}