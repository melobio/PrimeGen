import type { CompletedProtocolAnalysis } from '@xpcr/shared'

export function getLiquidDisplayName(
  analysis: CompletedProtocolAnalysis,
  liquidId: string
): CompletedProtocolAnalysis['liquids'][number]['displayName'] {
  const liquidDisplayName = (analysis?.liquids ?? []).find(
    liquid => liquid.id === liquidId
  )?.displayName
  return liquidDisplayName ?? ''
}
